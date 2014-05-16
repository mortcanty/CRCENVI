; docformat = 'rst'
; plr_reclass_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO plr_reclass_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'PLR Reclassify', $
      REF_VALUE = 'Overlay Classes', $
      EVENT_PRO = 'plr_reclass_run', $
      UVALUE = 'PLR_RECLASS',$
      POSITION = 'after'
END

;+
; :Description:
;       ENVI extension for postclassification with
;       Probabilistic Label Relaxation
;       Processes a rule image (class membership
;       probabilities), outputs a new classification 
;       file
; :Params:
;      event:  in, required 
;         if called from the ENVI menu                  
; :Uses:
;      ENVI::     
;      COYOTE
; :Author:
;      Mort Canty (2013)       
;-
pro plr_reclass_run, event

COMPILE_OPT IDL2

catch, theError
if theError ne 0 then begin
   print, 'Error reading file: '+fileName
   return
endif

print, '-----------------------------------------------'
print, 'Probabilistic Label Relaxation Reclassification'
print, systime(0)
print, '-----------------------------------------------'

lookup_flag = 0
envi_select, title='Choose classification file (Optional)', /no_spec, /no_dims, fid=fid, dims=dims, pos=pos
if (fid ne -1) then begin
   envi_file_query, fid, file_type=file_type
   if file_type ne 3 then begin
      error = dialog_message('Not a classification file',/error)
      return
   end
   envi_file_query, fid, lookup=lookup
end else begin
   void = dialog_message('Will use default color lookup table',/information)
   lookup_flag = 1
endelse

envi_select, title='Choose probabilities (rule) file', /no_spec, fid=fid, dims=dims, pos=pos
if (fid eq -1) then begin
   print,'canceled'
   return
endif
envi_file_query, fid, bnames=class_names, xstart=xstart, ystart=ystart
map_info = envi_get_map_info(fid=fid)
map_info.mc[0:1]= map_info.mc[0:1] - [dims[1],dims[3]]
n_cols = dims[2]-dims[1]+1
n_rows = dims[4]-dims[3]+1
n_classes = n_elements(pos)
prob_image = fltarr(n_cols,n_rows,n_classes)
for i=0,n_classes-1 do begin
    temp= envi_get_data(fid=fid,dims=dims,pos=pos[i])
    prob_image[*,*,i] = temp/255.0
endfor
query=dialog_message('Include dummy Unclassified class',/question)
if query eq 'Yes' then u_flag=0 else u_flag=1
base = widget_auto_base(title='Output PLR reclassification')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif
progressbar = Obj_New('cgprogressbar', /cancel,$
              title='Reclassifying...',xsize=250,ysize=20)
progressbar->start
class_image = bytarr(n_cols,n_rows)
for i=0L,n_cols-1 do begin
   if progressbar->CheckCancel() then begin
      print,'Reclassification aborted'
      progressbar->Destroy
      return
   endif
   pct=i*100/n_cols
   progressbar->Update,pct
   for j=0L,n_rows-1 do begin
       class_image[i,j]=(where(prob_image[i,j,*] eq max(prob_image[i,j,*])))[0]
       if (u_flag eq 0) and (max(prob_image[i,j,*]) eq 0) then class_image[i,j]=-1
   endfor
endfor
progressbar->destroy

; output
if u_flag eq 0 then begin
   class_names=['unclassified',class_names]
   n_classes=n_classes+1
   class_image=class_image+1B
endif
if lookup_flag then lookup = class_lookup_table(indgen(n_classes))

if (result.outf.in_memory eq 1) then begin
   envi_enter_data,class_image,file_type=3, $
                   map_info=map_info, $
                   bnames=['PLR Re-classification'],$
                   num_classes=n_classes, $
                   xstart=xstart, $
                   ystart=ystart, $
                   class_names=class_names, $
                   lookup=lookup
   print, 'PLR reclassification image written to memory'
endif else begin
   openw, unit, result.outf.name, /get_lun
   writeu, unit, class_image
   envi_setup_head ,fname=result.outf.name, ns=n_cols, $
                    nl=n_rows, nb=1, $
                    file_type=3,map_info=map_info,$
                    num_classes=n_classes, $
                    xstart=xstart, $
                    ystart=ystart, $
                    class_names=class_names, $
                    lookup=lookup, $
                    data_type=1, interleave=0, /write, $
                    bnames=['PLR Re-classification']
   print, 'PLR reclassification file created ', result.outf.name
   print, 'done'
   close, unit
   free_lun, unit
endelse

end