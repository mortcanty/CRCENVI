; docformat = 'rst'
; segment_class_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO segment_class_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Class Segmentation Image', $
      REF_VALUE = 'Buffer Zone Image', $
      EVENT_PRO = 'segment_class_run', $
      UVALUE = 'SEGMENT_CLASS',$
      POSITION = 'after'
END

;+
; :Description:
;       ENVI extension for segmentation of a classified 
;       image. Replaces ENVI built-in
; :Params:
;      event:  in, required 
;         if called from the ENVI menu                 
; :Uses:
;       ENVI::
;       SEGMENT_CLASS
; :Author:
;      Mort Canty (2013) 
pro segment_class_run, event

COMPILE_OPT IDL2

print, '--------------------------'
print, 'Classification Segmentation'
print, systime(0)
print, '--------------------------'

catch, theError
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
   return
endif

envi_select, title='Choose classification file', /no_spec, /no_dims, $
                    fid=fid, dims=dims, pos=pos
if (fid ne -1) then begin
   envi_file_query, fid, file_type=file_type,class_names=class_names,$
             data_type=data_type,fname=fname,xstart=xstart,ystart=ystart
   if file_type ne 3 then begin
      error = dialog_message('Not a classification file',/error)
      return
   endif       
end else begin
   print,'Cancelled'
   return
endelse 

num_cols = dims[2]-dims[1]+1
num_rows = dims[4]-dims[3]+1
 
climg = envi_get_data(fid=fid,dims=dims,pos=pos) 

; tie point
map_info = envi_get_map_info(fid=fid)
envi_convert_file_coordinates, fid, dims[1], dims[3], ee, nn, /to_map
map_info.mc[2:3]= [ee,nn]

base = widget_auto_base(title='Select classes')  
list = class_names[1:*] 
wm = widget_multi(base, list=list, uvalue='list', /auto)  
result = auto_wid_mng(base)  
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif  
; the K selected classes
cl = where(result.list eq 1,K)+1 

base = widget_auto_base(title='Minimum size')
wg = widget_sslider(base, title='SIze', min=0, max=500, $
  value=100, dt=2, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif
min_size = result.slide

base = widget_auto_base(title='Connectivity')  
list = ['    Four     ', '    Eight     ']  
sb = widget_base(base, /row)  
wt = widget_toggle(sb, uvalue='toggle', list=list, /auto)  
result = auto_wid_mng(base)  
if (result.accept eq 0) then all_neighbors=0 $  
else all_neighbors = result.toggle  

base = widget_auto_base(title='Segment Output')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif

; zero the edges
climg[0,*]=0
climg[num_cols-1,*]=0
climg[*,0]=0
climg[*,num_rows-1]=0     
        
; segment first class
limg = label_region(climg eq cl[0],/ULONG, $
        all_neighbors=all_neighbors)
; segment others, if any
for i = 1, K - 1 do begin
  mxlb   = max(limg)
  climgi = climg eq cl[i]
  limg   = temporary(limg) + $
        label_region(climgi, /ULONG, $
         all_neighbors=all_beighbors) $
        + climgi*mxlb
endfor
; re-label segments larger than min_size
indices = where(histogram(limg,reverse_indices=r) $
                           ge min_size,count)
limg=limg*0
if count gt 0 then for i=1L,count-1 do begin
   j = indices[i]       ;label of segment
   p = r[r[j]:r[j+1]-1] ;pixels in segment 
   limg[p]=i            ;new label
endfor

if (result.outf.in_memory eq 1) then $
   envi_enter_data, long(limg), $
          file_type=0, $
          map_info=map_info, $
          xstart=xstart+dims[1], $
          ystart=ystart+dims[3] $
else begin
   openw, unit, result.outf.name, /get_lun
   writeu, unit, long(limg)
   envi_setup_head ,fname=result.outf.name, ns=num_cols, nl=num_rows, nb=1, $
         data_type=3, $
         interleave=0, $
         file_type=0, $
         map_info=map_info, $
         xstart=xstart+dims[1], $
         ystart=ystart+dims[3], $
         descrip='Segmentation of classification' + result.outf.name, $
        /write
      print, 'File created ', result.outf.name
      free_lun, unit
endelse          

end