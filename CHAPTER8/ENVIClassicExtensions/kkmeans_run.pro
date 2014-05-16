; docformat = 'rst'
; kkmeans_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO kkmeans_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Kernel K-means (CUDA)', $
      REF_VALUE = 'K-Means', $
      EVENT_PRO = 'kkmeans_run', $
      UVALUE = 'KKMEANS',$
      POSITION = 'after'
END

;+
; :Description:
;       Performs kernel k-means clustering
;       using a Gaussian kernel::
;           Shawe-Taylor, J. and Cristianini, N. (2004).
;           Kernel Methods for Pattern Analysis. 
;           Cambridge University Press.
; :Params:
;       event:  in, optional 
;          required if called from ENVI                 
; :Uses:
;       ENVI::
;       GAUSSKERNEL_MATRIX::
;       GPUGAUSSKERNEL_MATRIX::  
;       CLASS_LOOKUP_TABLE::            
;       COYOTE::
;       GPULIB
; :Author:
;       Mort Canty (2013)      
;-
pro kkmeans_run, event

COMPILE_OPT IDL2

print, '-------------------'
print, 'Kernel K-means'
print, systime(0)
print, '-------------------'

;seed = 123L

catch, theError
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
   progressbar->destroy
   return
endif

envi_select, title='Choose multispectral image', $
             fid=fid, dims=dims,pos=pos
if (fid eq -1) then begin
   print, 'cancelled'
   return
endif
envi_file_query, fid, fname=fname,xstart=xstart,ystart=ystart
num_cols = dims[2]-dims[1]+1
num_rows = dims[4]-dims[3]+1
num_bands = n_elements(pos)
num_pixels = num_cols*num_rows
print, 'input file '+fname

; sample size
base = widget_auto_base(title='Sample size')
wg = widget_sslider(base, title='Samples', min=500, max=2000, $
  value=1000, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif 
m = long(result.slide)

; number of clusters
base = widget_auto_base(title='Number of Classes')
wg = widget_sslider(base, title='Classes', min=2, max=15, $
  value=6, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif
K = fix(result.slide)
print, 'Number of Classes', K

; nscale parameter for Gaussian kernel
base = widget_auto_base(title='NSCALE')  
we = widget_param(base, dt=4, field=3, floor=0.,xsize= 50,$  
  default=1.0, uvalue='param', /auto)  
result = auto_wid_mng(base)  
if (result.accept ne 0) then nscale = float(result.param) $
   else nscale = 1.0 

; export cluster to ROIs
query=dialog_message('Export clusters to ROIs',/question)
if query eq 'Yes' then roi_flag=1 else roi_flag=0

; output destination
base = widget_auto_base(title='K-means Output')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
  print, 'Output cancelled'
  return
endif

; map tie point
map_info = envi_get_map_info(fid=fid)
envi_convert_file_coordinates, fid, $
   dims[1], dims[3], e, n, /to_map
map_info.mc[2:3]= [e,n]

; image data matrix
GG = fltarr(num_bands,num_pixels)

progressbar = Obj_New('cgprogressbar',/cancel, $
              title='Kernel K-Means',xsize=300,ysize=20)
progressbar->start

; read in image 
for i=0,num_bands-1 do $
  GG[i,*] =  envi_get_data(fid=fid,dims=dims,pos=pos[i])

; training data matrix
indices = randomu(seed,m,/long) mod num_pixels
G = GG[*,indices]

; do not use zero data --------------------
tmp = total(abs(G),1)
idx = where(tmp gt 0, count)
if count gt 0 then begin
   G = G[*,idx]
   m = n_elements(idx)
endif    
; -----------------------------------------

; uncentered radial basis kernel matrix
if gpu_detect() then begin
   print,'running CUDA ...'
   G_gpu = gpuputarr(G)
   K_gpu = gpugausskernel_matrix(G_gpu,gma=gma,nscale=nscale)
   KK = gpugetarr(K_gpu)
   gpufree,K_gpu
end else begin
   print,'CUDA not available ...'
   KK = gausskernel_matrix(G,gma=gma,nscale=nscale)
endelse   

print,'GMA: ' + strtrim(gma,2)

; initial (random) class labels
labels = ceil(randomu(seed,m)*K)-1
; iteration
change = 1
iter = 0
ones = fltarr(1,m)+1
while change and (iter lt 100) do begin
   progressbar->Update,iter 
   change = 0
   U = fltarr(m,K)
   for i=0,m-1 do U[i,labels[i]] = 1
   M1 = diag_matrix(1/(total(U,1)+1))
   MU = M1##U
; Z is a K by m array 
   Z = ones##diag_matrix(MU##KK##transpose(MU)) $
                   - 2*KK##transpose(MU)
   _ = min(Z,labels1,dimension=1)
   labels1 = labels1 mod K
   if total(labels1 ne labels) then change=1
   labels=labels1
   if progressbar->CheckCancel() then begin
      print,'interrupted...'
      iter=99
   endif
   iter++
endwhile
print,'Iterations: '+strtrim(iter,2)
if iter eq 100 then print,'Interrupted or no convergence after 100 iterations'
progressbar->destroy

; export clusters to ROIs 
labels = labels+1
c_names = strarr(K+1)
c_names[0] = 'unclassified'
for i=1,K do c_names[i]='cluster: '+strtrim(i,2)
if roi_flag then begin
   envi_delete_rois,envi_get_roi_ids()
   for i=1,K do begin
      roi_id = envi_create_roi(color=i,name=c_names[i],ns=num_cols,nl=num_rows)
      ind = where(labels eq i, count)
      if count gt 0 then begin
         xpts = indices[ind] mod num_cols
         ypts = indices[ind]/num_cols   
         envi_define_roi, roi_id, /point, xpts=xpts+dims[1], ypts=ypts+dims[3]
      endif
   endfor
endif

; classify image
start_time=systime(2)
progressbar = Obj_New('cgprogressbar',/cancel, $
           title='Classifiying',xsize=300,ysize=20)
progressbar->start 
if gpu_detect() then begin
; GPU variables
; A_gpu is K x num_cols
   gpuPutArr,(fltarr(num_cols)+1)##diag_matrix(MU##KK##transpose(MU)),A_gpu
   gpuPutArr,MU,MU_gpu
   gpuPutArr,GG,GG_gpu    ;data matrix
   _ = temporary(Z)
; loop over rows
   i=0L
   class_image = bytarr(num_cols,num_rows)             
   while i lt num_rows do begin
      pct=i*100/num_rows
      progressbar->Update,fix(pct)
;    ggi = GG[*,i*num_cols:(i+1)*num_cols-1]
      gpuView,GG_gpu,i*num_cols*num_bands,num_bands*num_cols,GGi_gpu
      gpuReform,GGi_gpu,num_bands,num_cols 
;    KKK = gausskernel_matrix(G,GGi,gma=gma) ; num_cols x m 
      KKK_gpu = gpugausskernel_matrix(G_gpu,GGi_gpu,gma=gma)      
;    Z is K x num_cols
;    Z = A - 2*transpose(KKK)##transpose(MU)
      Z_gpu = gpuMatrix_Multiply(MU_gpu,KKK_gpu,/atranspose,/btranspose)
      Z_gpu = gpuAdd(1.0,A_gpu,-2.0,Z_gpu,0.0,lhs=Z_gpu)
      Z = gpugetarr(Z_gpu)
      _ = min(Z,labels,dimension=1)
      gpuFree,[Z_gpu,KKK_gpu]
      class_image[*,i] = byte(labels mod K)+1
      if progressbar->CheckCancel() then begin
         print,'interrupted...'
         i=num_rows-1
      endif
      i++ ; next row
   endwhile
; tidy up
   gpuFree,[A_gpu,MU_gpu,GG_gpu,G_gpu]
end else begin
   i=0L
; A is (K by num_cols) array
   A = (fltarr(num_cols)+1)##diag_matrix(MU##KK##transpose(MU))
   class_image = bytarr(num_cols,num_rows)
   while i lt num_rows do begin
      pct=i*100/num_rows
      progressbar->Update,fix(pct)
;    get the ith row of the image
      GGi = GG[*,i*num_cols:(i+1)*num_cols-1]
;    (num_cols by m) array
      KKK = gausskernel_matrix(G,GGi,gma=gma)
;    Z is a (K by num_cols) array 
      Z = A - 2*transpose(KKK)##transpose(MU)
;    classify
      _ = min(Z,labels,dimension=1)
      class_image[*,i] = byte(labels mod K)+1
      if progressbar->CheckCancel() then begin
         print,'interrupted...'
         i=num_rows-1
      endif
      i++ ; next row
   endwhile
endelse    
progressbar->destroy
print,'elapsed time: '+ strtrim(systime(2)-start_time,2)+' seconds'

; write result to memory or disk
if (result.outf.in_memory eq 1) then begin
   envi_enter_data, class_image, $
                    file_type=3, $
                    map_info=map_info, $
                    bnames=['KKM('+fname+')'],$
                    num_classes=K+1, $
                    class_names=c_names, $
                    xstart=xstart+dims[1], $
                    ystart=ystart+dims[3], $
                    lookup=class_lookup_table(indgen(K+1))
   print, 'Result written to memory'
endif else begin
   openw, unit, result.outf.name, /get_lun
   writeu, unit, class_image
   envi_setup_head ,fname=result.outf.name, ns=num_cols, nl=num_rows, nb=1, $
       data_type=1, $
       interleave=2, $
       file_type=3, $
       map_info=map_info, $
       bnames=['FKM('+fname+')'],$
       num_classes=K+1, $
       class_names=c_names, $
       xstart=xstart+dims[1], $
       ystart=ystart+dims[3], $
       lookup=class_lookup_table(indgen(K+1)), $
       descrip='KKM clustering of ' + result.outf.name, $
       /write
   print, 'File created ', result.outf.name
   free_lun, unit
endelse

end
