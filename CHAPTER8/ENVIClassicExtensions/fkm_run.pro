; docformat = 'rst'
; fkm_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO fkm_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Fuzzy K-Means', $
      REF_VALUE = 'K-Means', $
      EVENT_PRO = 'fkm_run', $
      UVALUE = 'FKMM',$
      POSITION = 'after'
END

;+
; :Description:
;       ENVI extension for fuzzy K-means clustering
;       with sampled data  
; :Params:
;       event:  in, optional 
;          required if called from ENVI                 
; :Uses:
;       ENVI::
;       FMK::  
;       CLUSTER_FKM (Modified distance clusterer 
;       from IDL library)::
;       CLASS_LOOKUP_TABLE            
; :Author:
;       Mort Canty (2013)      
;-
pro fkm_run, event

COMPILE_OPT IDL2

print, '---------------------------------'
print, 'Fuzzy-K-Means Clustering'
print, systime(0)
print, '---------------------------------'

catch, theError
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
   return
endif

envi_select, title='Choose multispectral image for clustering', $
             fid=fid, dims=dims,pos=pos
if (fid eq -1) then begin
   print, 'cancelled'
   return
endif

envi_file_query, fid, fname=fname, xstart=xstart, ystart=ystart
print, 'Selected image: ',fname
print, 'Selected bands: ',pos + 1

; tie point
map_info = envi_get_map_info(fid=fid)
envi_convert_file_coordinates, fid, dims[1], dims[3], e, n, /to_map
map_info.mc[2:3]= [e,n]

; number of clusters
base = widget_auto_base(title='Number of Classes')
wg = widget_sslider(base, title='Classes', min=2, max=15, $
  value=6, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif
K = byte(result.slide)
print, 'Number of clusters', K

query=dialog_message('Export clusters to ROIs',/question)
if query eq 'Yes' then roi_flag=1 else roi_flag=0

; output destination
base = widget_auto_base(title='FKM Output')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
  print, 'cancelled'
  return
endif

c_names = strarr(K+1)
c_names[0] = 'unclassified'
for i=1,K do c_names[i]='cluster: '+strtrim(i,2)

num_cols = dims[2]-dims[1]+1L
num_rows = dims[4]-dims[3]+1L
num_pixels = (num_cols*num_rows)
num_bands = n_elements(pos)

widget_control, /hourglass

if num_pixels gt 100000 then $
; random sample of 100000 pixels
   indices = randomu(seed, 100000, /long) mod num_pixels $
; all pixels
else indices = lindgen(num_pixels)

G = fltarr(num_bands,n_elements(indices))
image = fltarr(num_bands,num_pixels)
for i=0,num_bands-1 do begin
    temp= envi_get_data(fid=fid,dims=dims,pos=pos[i])
    image[i,*] = temp
    G[I,*]=temp[indices]
endfor

; initialize memberships at random
U = randomu(seed,n_elements(indices),K)
for j=0,K-1 do U[*,j]=U[*,j]/total(U,2)

; run fuzzy-K-means
print, 'FKM-algorithm running on ',n_elements(indices),' pixels ...'
FKM, G, U, Ms

; export clusters to ROIs
if roi_flag then begin
   envi_delete_rois, envi_get_roi_ids()
   labels = cluster_fkm(G, Ms) + 1B
   if n_elements(lables) eq 1 then return
   for i=1,K do begin
      roi_id = envi_create_roi(color=i+1, name=c_names[i], ns=num_cols, nl=num_rows)
      ind = where(labels eq i,count)
      if count gt 0 then begin
         xpts = indices[ind] mod num_cols
         ypts = indices[ind]/num_cols
      endif
      envi_define_roi, roi_id, /point, xpts=xpts, ypts=ypts
    endfor
endif

; classify the image
print, 'classifying...'
labels = cluster_fkm(image, Ms) + 1B
if n_elements(labels) eq 1 then return

; write result to memory or disk
if (result.outf.in_memory eq 1) then begin
   envi_enter_data, byte(reform(labels,num_cols,num_rows)), $
                    file_type=3, $
                    map_info=map_info, $
                    bnames=['FKM('+fname+')'],$
                    num_classes=K+1, $
                    class_names=c_names, $
                    xstart=xstart+dims[1], $
                    ystart=ystart+dims[3], $
                    lookup=class_lookup_table(indgen(K+1))
   print, 'Result written to memory'
endif else begin
   openw, unit, result.outf.name, /get_lun
   writeu, unit, byte(reform(labels,num_cols,num_rows))
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
       descrip='FKM clustering of ' + result.outf.name, $
       /write
   print, 'File created ', result.outf.name
   free_lun, unit
endelse

end



