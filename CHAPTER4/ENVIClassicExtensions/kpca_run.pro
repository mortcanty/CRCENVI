; docformat = 'rst'
; kpca_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO kpca_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Kernel PCA (CUDA)', $
      REF_VALUE = 'Inverse PC Rotation', $
      EVENT_PRO = 'kpca_run', $
      UVALUE = 'KPCA',$
      POSITION = 'after' 
END

;+
; :Description:
;      Performs kernel principal components analysis
;      using a Gaussian kernel::      
;         Nielsen, A. A. and Canty, M. J. (2008). 
;         Kernel principal component analysis
;         for change detections. In SPIE Europe 
;         Remote Sensing Conference, Cardiff,
;         Great Britain, 15-18 September, volume 7109.       
; :Params:
;      event:  in, required 
;         if called from the ENVI menu                  
; :Uses:
;      ENVI, CENTER, KERNEL_MATRIX, GPUKERNEL_MATRIX, GPUKERNELPROJECT, COYOTE, GPULIB
; :Author:
;      Mort Canty (2013)       
;-
pro KPCA_run, event 

COMPILE_OPT IDL2

seed = 12345L
  
print, '-------------------'
print, 'Kernel PCA'
print, systime(0)
print, '-------------------'

catch, theError
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
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
base = widget_auto_base(title='Sample size, 0 for k-means')
wg = widget_sslider(base, title='Samples', min=0, max=5000, $
  value=2000, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif 
m = long(result.slide)
if m gt 0 then $
   print,'Training on '+strtrim(m,2)+' samples' $
else print, 'Using k-means'    

; number of PCs
base = widget_auto_base(title='Numer of PCs')
wg = widget_sslider(base, title='PCs', min=5, max=50, $
  value=10, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif 
num_pcs = long(result.slide)

; nscale parameter for Gaussian kernel
base = widget_auto_base(title='NSCALE')  
we = widget_param(base, dt=4, field=3, floor=0.,xsize= 50,$  
  default=1.0, uvalue='param', /auto)  
result = auto_wid_mng(base)  
if (result.accept ne 0) then nscale = float(result.param) $
   else nscale = 1.0 
   
query=dialog_message('Center on training means',/question)
if query eq 'Yes' then ct=1 else ct=0   

; output destination
base = widget_auto_base(title='KPCA Output')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result1 = auto_wid_mng(base)
if (result1.accept eq 0) then begin
  print, 'Output cancelled'
  return
end

progressbar = Obj_New('cgprogressbar', Text='Initializing...',$
             title='Kernel PCA ',xsize=300,ysize=20,/nocancel)            
progressbar->start

start_time = systime(2) 

; image data matrix
GG = fltarr(num_bands,num_pixels)
for i=0,num_bands-1 do $
 GG[i,*] = envi_get_data(fid=fid,dims=dims,pos=pos[i])

if m gt 0 then begin
; random training data matrix  
   indices = randomu(seed,m,/long) mod num_pixels
   G = GG[*,indices]
endif else begin
; k-means training matrix
   m = 100
   G = CLUST_WTS(GG,N_CLUSTERS=m)
   print,'Elapsed time for k-means: '+strtrim(systime(2)-start_time,2)+' seconds'   
endelse  

 start_time = systime(2)  

; centered radial basis kernel matrix
if gpu_detect() then begin
   G_gpu = gpuputarr(G)
   K_gpu = gpukernel_matrix(G_gpu,gma=gma,nscale=nscale)
   K = gpugetarr(K_gpu)
   gpufree,[G_gpu,K_gpu]
end else begin
   K = kernel_matrix(G,gma=gma,nscale=nscale)
endelse     
K = center(K)

print,'GMA: ' + strtrim(gma,2)

progressbar->Update,0

; eigenvalues and eigenvectors of centered kernel matrix
print, 'diagonalizing kernel...'
lambda = la_eigenql(K,/double,eigenvectors=V,range=[m-num_pcs,m-1])
idx = reverse(sort(lambda))
lambda = lambda[idx]
V = V[*,idx]

print,'First '+strtrim(num_pcs,2)+' eigenvalues of the kernel matrix'
print, lambda
print, 'Fraction of variance explained in first 3 PCS'
print, total(lambda[0:2])/total(lambda)
   
envi_plot_data,findgen(num_pcs)+1,lambda[0:num_pcs-1],$
   plot_title='Kernel PC File: '+file_basename(fname),$
   title='Kernel PCA',$
   ytitle='Eigenvalue',$
   xoff=100,yoff=100

print,'Elapsed time for training: '+strtrim(systime(2)-start_time,2)+' seconds'
progressbar->destroy

;if num_bands eq 2 then begin
;   query=dialog_message('Generate contour image',/question)
;   if query eq 'Yes' then begin
;;    contours of first num_pcs principal components
;      FF = fltarr(2,400,400)
;      FF[1,*,*] = transpose(lindgen(400))##(fltarr(400)+1)/399
;      FF[0,*,*] = transpose(FF[1,*,*])
;      FF = reform(FF,2,160000L)
;      for i=0,1 do $
;         FF[i,*] = min(G[i,*]) + FF[i,*]*(max(G[i,*])-min(G[i,*]))
;      alpha = float(diag_matrix(1/sqrt(lambda))##V)   
;      if gpuKernelProject(alpha, gma, G, FF, 400, 400, image) then begin
;         envi_enter_data, reverse(image,2)      
;;       export training data to ROI
;         xpts = round( (G[0,*]-min(G[0,*]))*400/(max(G[0,*])-min(G[0,*])) )
;         ypts = round( (G[1,*]-min(G[1,*]))*400/(max(G[1,*])-min(G[1,*])) ) 
;         envi_delete_rois,envi_get_roi_ids()
;         roi_id = envi_create_roi(color=0,name='training data',ns=400,nl=400)
;         envi_define_roi, roi_id, /point, xpts=xpts[*], ypts=400-ypts[*] 
;      endif   
;   endif           
;endif

start_time = systime(2) 

; dual variables (normalized eigenvectors)
alpha = float(diag_matrix(1/sqrt(lambda))##V) 
   
if not gpuKernelProject(alpha,  gma, G, GG, num_cols, num_rows, image, center_train=ct) then $
       message, 'projection aborted'

print,'Elapsed time for projection: '+strtrim(systime(2)-start_time,2)+' seconds'
         
; map tie point
map_info = envi_get_map_info(fid=fid)
envi_convert_file_coordinates,fid,dims[1],dims[3],e,n,/to_map
map_info.mc = [0D,0D,e,n]

; write to memory or file
bnames='kernel PCA '+strtrim(lindgen(num_pcs)+1,2)
if (result1.outf.in_memory eq 1) then begin
   envi_enter_data, image, $
      map_info=map_info, $
      bnames=bnames, $
      xstart=xstart+dims[1], ystart=ystart+dims[3], $
      descrip='kernel PCA: '+file_basename(fname)
   print, 'PCs written to memory'
end else begin
   openw, unit, result1.outf.name, /get_lun
   writeu, unit, image
   envi_setup_head,fname=result1.outf.name, ns=num_cols, nl=num_rows, nb=num_pcs, $
                    data_type=4, $
                    interleave=0, $
                    file_type=0, $
                    map_info=map_info, $
                    xstart=xstart+dims[1], $
                    ystart=ystart+dims[3], $
                    bnames=bnames,$            
                    descrip='kernel PCA: '+file_basename(fname), $
                    /write,/open
   print, 'File created ', result1.outf.name
   free_lun, unit
endelse

end
