; docformat = 'rst'
; plr_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO plr_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'PLR', $
      REF_VALUE = 'Overlay Classes', $
      EVENT_PRO = 'PLR_run', $
      UVALUE = 'PLR',$
      POSITION = 'after', $
      /SEPARATOR
END

;+
; :Description:
;       ENVI extension for postclassification with
;       Probabilistic Label Relaxation.
;       Processes a rule image (class membership
;       probabilities), outputs a new rule image::  
;          Richards, J. A. and Jia, X. (2006). Remote 
;          Sensing Digital Image Analysis ;(4th Ed.). 
;          Springer.
; :Params:
;      event:  in, required 
;         if called from the ENVI menu                  
; :Uses:
;      ENVI::     
;      COYOTE
; :Author:
;      Mort Canty (2013)       
;-
pro plr_run, event

COMPILE_OPT IDL2

catch, theError
if theError ne 0 then begin
   print, 'Error reading file: '+fileName
   return
endif

print, '---------------------------------'
print, 'Probabilistic Label Relaxation'
print, systime(0)
print, '---------------------------------'

envi_select, title='Choose probabilities image', /no_spec, fid=fid, dims=dims, pos=pos
if (fid eq -1) then begin
   print,'canceled'
   return
endif
envi_file_query, fid, bnames=bnames,xstart=xstart, ystart=ystart, data_type=data_type
map_info = envi_get_map_info(fid=fid)
map_info.mc[0:1]= map_info.mc[0:1] - [dims[1],dims[3]]
n_cols = dims[2]-dims[1]+1L
n_rows = dims[4]-dims[3]+1L
n_classes = n_elements(pos)
prob_image = fltarr(n_cols,n_rows,n_classes)
for i=0,n_classes-1 do begin
    temp= envi_get_data(fid=fid,dims=dims,pos=pos[i])
    if data_type eq 1 then prob_image[*,*,i] = temp/255.0 $
        else prob_image[*,*,i] = temp
endfor
base = widget_auto_base(title='Number of Iterations')
wg = widget_sslider(base, title='Iterations', min=1, max=10, $
  value=3, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print,'cancelled'
   return
endif
nIter = fix(result.slide)
base = widget_auto_base(title='PLR Output')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
  print, 'cancelled'
  return
endif
; compatibility matrix
widget_control,/hourglass
n_samples=(n_cols-1)*(n_rows-1)
samplem = reform(prob_image[0:n_cols-2,0:n_rows-2,*],n_samples,n_classes)
samplen = reform(prob_image[1:n_cols-1,0:n_rows-2,*],n_samples,n_classes)
sampleu = reform(prob_image[0:n_cols-2,1:n_rows-1,*],n_samples,n_classes)
max_samplem = max(samplem,dimension=2)
max_samplen = max(samplen,dimension=2)
max_sampleu = max(sampleu,dimension=2)
Pmn = fltarr(n_classes, n_classes)
progressbar = Obj_New('cgprogressbar', /cancel,$
              title='PLR: compatibility matrix...',xsize=250,ysize=20)
progressbar->start
for j=0L, (n_samples-1) do begin
   if j mod 1000 eq 0 then begin
      if progressbar->CheckCancel() then begin
         print,'PLR aborted'
         progressbar->Destroy
         return
      endif
      pct=j*100/n_samples
      progressbar->Update,pct
   endif
   samplem1 = (where(samplem[j, *] eq max_samplem[j],countm))[0]
   samplen1 = (where(samplen[j, *] eq max_samplen[j],countn))[0]
   if countm gt 0 and countn gt 0 then Pmn[samplen1, samplem1] = Pmn[samplen1, samplem1] + 1
   sampleu1 = (where(sampleu[j, *] eq max_sampleu[j],countu))[0]
   if countm gt 0 and countu gt 0 then Pmn[sampleu1, samplem1] = Pmn[sampleu1, samplem1] + 1
endfor
progressbar->Destroy
for j=0, n_classes-1 do begin
   n = total(Pmn[*,j])
   if (n ne 0) then Pmn[*,j] = (Pmn[*,j]) / n
endfor
print,'Compatibility matrix'
print,Pmn
; start iterations
temp = prob_image*0.0
for iter=1,nIter do begin
; relaxation
   Pm = fltarr(n_classes)
   Pn = fltarr(n_classes)
   progressbar = Obj_New('cgprogressbar', /cancel,$
         title='PLR: pass '+strtrim(iter,2)+' of '+strtrim(nIter,2)+' ...',xsize=250,ysize=20)
   progressbar->start
   for i=1L, n_cols-2 do begin
      if progressbar->CheckCancel() then begin
         print,'PLR aborted'
         progressbar->Destroy
         return
      endif
      pct = i*100/n_cols
      progressbar->Update,pct
      for j=1L, n_rows-2 do begin
         Pm[*] = prob_image[i, j,*]
         Pn[*] = prob_image[i-1, j, *]/4
         Pn[*] = Pn[*] + prob_image[i+1, j, *]/4
         Pn[*] = Pn[*] + prob_image[i, j-1, *]/4
         Pn[*] = Pn[*] + prob_image[i, j+1, *]/4
         if max(Pm) eq 0 then Pm_new = Pm  $
            else Pm_new = Pm*(Pmn##transpose(Pn))/(Pm##Pmn##transpose(Pn))[0]
         temp[i, j,*] = Pm_new
      endfor
   endfor
   prob_image = temp
   progressbar->destroy
endfor
; end iterations, save results
if (result.outf.in_memory eq 1) then begin
   envi_enter_data, bytscl(prob_image,min=0.0,max=1.0),$
     map_info=map_info, $
     xstart=xstart, $
     ystart=ystart, $
     bnames=bnames
   print, 'Result written to memory'
endif else begin
   openw, unit, result.outf.name, /get_lun
   writeu, unit, bytscl(prob_image,min=0.0,max=1.0)
   free_lun, unit
   envi_setup_head, fname=result.outf.name, ns=n_cols, $
       nl=n_rows, nb=n_classes, $
       data_type=1, interleave=0, /write, $
       map_info=map_info, $
       bnames=bnames, $
       xstart=xstart, $
       ystart=ystart, $
       descrip='PLR probabilities'
   print, 'File created ', result.outf.name
   close, unit
endelse

end