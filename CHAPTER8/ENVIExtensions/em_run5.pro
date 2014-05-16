;docformat = 'rst'
; em_run5.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

pro em_run5_extensions_init
   compile_opt idl2
   e = envi(/current)
   e.addextension, 'Gaussian Mixture Clustering (CUDA)', 'em_run5', path='Chapter8'
end

function cov_to_corr, A
   n = n_elements(A[*,0])
   v = diag_matrix(1/sqrt(diag_matrix(A)))##(fltarr(n,n)+1)
   return, transpose(v)*A*v
end

;+
; :Description:
;       ENVI extension for Gaussian mixture clustering 
;       at different scales using DWT compression,
;       simulated annealing, Gauss-Markov
;       random fields and the EM algorithm
; :Params:
;       NONE               
; :Uses:
;       ENVI::
;       EM::
;       GPUEM::  
;       DWT__DEFINE:: 
;       CLASS_LOOKUP_TABLE::           
; :Author:
;       Mort Canty (2013)      
;-
pro em_run5

COMPILE_OPT IDL2
common sd, seed

print, '---------------------------------'
print, 'Scaled EM Clustering'
print, systime(0)
print, '---------------------------------'

;seed = 123L

catch, err
   if err ne 0 then begin
      catch, /cancel
      if e5 ne !null then $
         e5.reporterror, 'Error: ' + !error_state.msg $
      else $
         message, !error_state.msg, /continue, /noname
      message, /reset
      return
   endif
   
e5 = envi(/current)
   if e5 eq !null then $
      message, 'This extension requires an interactive ENVI session 

envi_select, title='Choose multispectral image', $
             fid=fid, dims=dims,pos=pos
if (fid eq -1) then begin
   print, 'cancelled'
   return
endif

envi_file_query, fid, fname=fname, xstart=xstart, ystart=ystart, ns=ns, nl=nl
print, 'Selected image: ',fname
print, 'Selected band(s): ',strtrim(pos + 1,2)

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
print, 'Number of Classes: ', strtrim(K,2)

; initial scale
base = widget_auto_base(title='Maximum scaling factor')
wg = widget_sslider(base, title=' ', min=0, max=3, $
  value=2, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then return
max_scale = fix(result.slide)
print, 'Maximum scale: ', strtrim(max_scale,2)

base = widget_auto_base(title='Minimum scaling factor')
wg = widget_sslider(base, title=' ', min=0, max=3, $
  value=0, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then return
min_scale = fix(result.slide)
min_scale = min_scale < max_scale
print, 'Minimum scale: ', strtrim(min_scale,2)

; tie point
map_info = envi_get_map_info(fid=fid)
envi_convert_file_coordinates, fid, dims[1], dims[3], e, n, /to_map
map_info.mc[2:3]= [e,n]

base = widget_auto_base(title='Clustering parameters')
list = ['T0', 'Beta']
vals = [0.5,0.5]
we = widget_edit(base,  list=list, uvalue='edit', $
      vals=vals, field= 2, dt=4, /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print,'cancelled'
   return
endif
T0 = result.edit[0]
beta = result.edit[1]

query=dialog_message('Export clusters to ROIs',/question)
if query eq 'Yes' then roi_flag=1 else roi_flag=0

; output destination
base = widget_auto_base(title='EM Output')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result1 = auto_wid_mng(base)
if (result1.accept eq 0) then begin
  print, 'Output cancelled'
  return
endif

base = widget_auto_base(title='Output class membership probs')
sb = widget_base(base, /row, /frame)
wp = widget_outf(sb, uvalue='outf', /auto)
result2 = auto_wid_mng(base)
if (result2.accept eq 0) then begin
  print, 'Output of class membership probs cancelled'
  prob_flag=0
end else prob_flag=1

c_names=strarr(K+1)
fhvs = fltarr(K)
c_names[0]='unclassified'
for i=1,K do c_names[i]='cluster: '+strtrim(i,2)

num_cols = dims[2]-dims[1]+1L
num_rows = dims[4]-dims[3]+1L
num_pixels = (num_cols*num_rows)
num_bands = n_elements(pos)

widget_control, /hourglass

; filter
DWTbands = objarr(num_bands)
for i=0,num_bands-1 do begin
   print,'DWT of band '+strtrim(i+1,2)+'...'
   DWTbands[i] = Obj_New('DWT',envi_get_data(fid=fid,dims=dims,pos=pos[i]))
   for c=1,max_scale do DWTbands[i]->filter
endfor

s_num_cols = DWTBands[0]->get_num_cols()
s_num_rows = DWTBands[0]->get_num_rows()
s_num_pixels = s_num_cols*s_num_rows
s_image = fltarr(num_bands,s_num_pixels)
for i=0,num_bands-1 do s_image[i,*] = DWTbands[i]->get_quadrant(0)

; initialize membership matrix 
U = randomu(seed,s_num_pixels,K)
den = total(U,2)
for j=0,K-1 do U[*,j]=U[*,j]/den

; initial clustering of maximally compressed image
window,11,xsize=600,ysize=400,title='Log Likelihood'
if gpu_detect() and (num_bands gt 1) then $
  gpuEM, s_image, U, Ms, Ps, Cs, T0=T0, wnd=11, $
      num_cols=s_num_cols, num_rows=s_num_rows, $
      beta=beta, pdens=pdens,status=status,/verbose $
else EM,  s_image, U, Ms, Ps, Cs, T0=T0, wnd=11, $
      num_cols=s_num_cols, num_rows=s_num_rows, $
      beta=beta, pdens=pdens,status=status,/verbose      
if status then message, 'EM algoritm failed'

; sort clusters wrt partition density
U = U[*,reverse(sort(pdens))]

; clustering at increasing scales
for scale=1,max_scale-min_scale do begin
; upsample class memberships (bilinear interpolation)
   U = reform(U,s_num_cols,s_num_rows,K,/overwrite)
   s_num_cols = s_num_cols*2
   s_num_rows = s_num_rows*2
   s_num_pixels = s_num_cols*s_num_rows
   U = shift(rebin(U,s_num_cols,s_num_rows,K),1,1,0)
   U = reform(U,s_num_pixels,K,/overwrite)
   
; expand the image
   s_image = fltarr(num_bands,s_num_pixels)
   for i=0,num_bands-1 do begin
      DWTBands[i]->invert
      s_image[i,*] = DWTbands[i]->get_quadrant(0)
   endfor
   
; cluster
   unfrozen = where(max(U,dimension=2) le 0.90, count)
   if count gt 0 then $
     if gpu_detect() and (num_bands gt 1) then $
       gpuEM,s_image,U,Ms,Ps,Cs,miniter=5,wnd=11, $
          unfrozen=unfrozen,T0=0.0, $
          num_cols=s_num_cols, num_rows=s_num_rows, $
          beta=beta, status=status,/verbose $
     else EM,s_image,U,M,Ps,Cs,miniter=5,wnd=11, $
            unfrozen=unfrozen, T0=0.0, $
            num_cols=s_num_cols,num_rows=s_num_rows, $
            beta=beta, status=status,/verbose  
   if status then return
endfor

if num_bands gt 1 then begin
   print, 'Cluster mean vectors'
   print, Ms[*,0:K-1]
   print, 'Cluster covariance matrices'
   for i=0, K-1 do begin
      print, strtrim(i+1,2)
      print, Cs[*,*,i] 
   endfor
   print, 'Cluster correlation matrices'
   for i=0, K-1 do begin
      print, strtrim(i+1,2)
      print, cov_to_corr( Cs[*,*,i] )
   endfor
end else begin
   print, 'Cluster means'
   print, transpose(Ms)
   print, 'Cluster variances'
   print, Cs
endelse  

; if image is still downscaled, upsample the class memberships
if min_scale gt 0 then  begin
   U = reform(U,s_num_cols,s_num_rows,K,/overwrite)
   s_num_cols=(2^min_scale)*s_num_cols
   s_num_rows=(2^min_scale)*s_num_rows
   s_num_pixels=s_num_cols*s_num_rows
   U = shift(rebin(U,s_num_cols,s_num_rows,K),2^min_scale,2^min_scale,0)
   U = reform(U,s_num_pixels,K,/overwrite)
endif

; classify the final expanded image
maxima = max(U,labels,dimension=2)
if n_elements(labels) eq 1 then goto, done
labels = byte((labels/s_num_pixels)+1)
; get the unclassified pixels ------;
indices = where(maxima eq 0,count)  ; outliers for which U = 0
if count gt 0 then labels[indices]=0; are set as unclassified
; ----------------------------------;
out_array = bytarr(num_cols,num_rows)
out_array[0:s_num_cols-1,0:s_num_rows-1] = reform(labels,s_num_cols,s_num_rows)

; export clusters to ROIs (100000 samples in all)
if roi_flag then begin
   indices = lindgen(s_num_pixels)
   envi_delete_rois,envi_get_roi_ids()
   for i=1,K do begin
      roi_id = envi_create_roi(color=i,name=c_names[i],ns=ns,nl=nl)
      ind = where(labels eq i, count)
      numsamples = fix((100000.0 * count)/s_num_pixels,type=3)
      if count gt 0 then begin
         sample = randomu(seed,numsamples,/long) mod count
         ind = ind[sample]
         xpts = indices[ind] mod s_num_cols
         ypts = indices[ind]/s_num_cols
         envi_define_roi, roi_id, /point, xpts=xpts+dims[1], ypts=ypts+dims[3]
      endif      
   endfor
endif

; output classification image to memory or file
if (result1.outf.in_memory eq 1) then begin
   envi_enter_data, out_array+0, $
                    file_type=3, $
                    map_info=map_info, $
                    xstart=xstart+dims[1], $
                    ystart=ystart+dims[3], $
                    bnames=['EM('+fname+')'],$
                    num_classes=K+1, $
                    class_names=c_names, $
                    lookup=class_lookup_table(indgen(K+1))
   print, 'Classification written to memory'
end else begin
   openw, unit, result1.outf.name, /get_lun
   writeu, unit, out_array
   envi_setup_head ,fname=result1.outf.name, ns=num_cols, nl=num_rows, nb=1, $
                    data_type=1, $
                    interleave=2, $
                    file_type=3, $
                    map_info=map_info, $
                    xstart=xstart+dims[1], $
                    ystart=ystart+dims[3], $
                    bnames=['EM('+fname+')'],$
                    num_classes=K+1, $
                    class_names=c_names, $
                    lookup=class_lookup_table(indgen(K+1)), $
                    descrip='scaled EM clustering of ' + fname, $
                    /write
   print, 'File created ', result1.outf.name
   free_lun, unit
endelse

if not prob_flag then goto, done

; Output membership probabilities to file

; convert to byte to save memory
U = bytscl(U,min=0.0,max=1.0)

openw, unit, result2.outf, /get_lun
for i=0,K-1 do begin
   out_array[0:s_num_cols-1,0:s_num_rows-1,*]=reform(U[*,i],s_num_cols,s_num_rows)
   writeu, unit, out_array
endfor
free_lun, unit
envi_setup_head, fname=result2.outf, ns=num_cols, nl=num_rows, nb=K, $
    data_type=1, $
    interleave=0, $
    file_type=0, $
    map_info=map_info, $
    xstart=xstart+dims[1], $
    ystart=ystart+dims[3], $
    bnames='probs: '+c_names, $
    descrip='EM probabilities for ' + fname, $
    /write
print, 'File created ', result2.outf

done:
for i=0,num_bands-1 do Obj_Destroy, DWTbands[i]

end