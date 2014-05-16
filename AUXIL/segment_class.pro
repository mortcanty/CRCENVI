; docformat = 'rst'
; segment_class.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details
;+
; :Description:
;      perform blobbing on a segmented or classified image
; :Params:
;      climg: in, required
;         classification image to be segmented
;      cl: in, required, type=array of integer            
;         classes to be included in segmentation
; :Keywords:
;      minsize: in, optional
;         minimum segment size (default 0)
;      all_neighbours: in, optional
;         set for 8-neighborhood
;         (default4-neighborhood)
; :Uses:
;      COYOTE         
; :Author:
;      Mort Canty (2009) 
;-
function segment_class, climg, cl, min_size=min_size, all_neighbors=all_neighbors

   if n_elements(all_neighbors) eq 0 then all_neighbors=0
   if n_elements(min_size) eq 0 then min_size=0
   
   num_cols = (size(climg))[1]
   num_rows = (size(climg))[2]
   K = n_elements(cl)

; zero the edges
   climg[0,*]=0
   climg[num_cols-1,*]=0
   climg[*,0]=0
   climg[*,num_rows-1]=0     
        
; segment first class
   limg = LABEL_REGION(climg EQ cl[0],/ULONG,all_neighbors=all_neighbors)
; segment others, if any
   progressbar = Obj_New('progressbar', Color='blue', Text=' ',$
                 title='Blobbing...',xsize=300,ysize=20)
   progressbar->start
   FOR i = 1, K - 1 DO BEGIN
     if progressbar->CheckCancel() then begin
        progressbar->destroy
        return,-1
     endif
     pct = (i)*100.0/K
     progressbar->Update,fix(pct),text=strtrim(i)  
     mxlb   = MAX(limg)
     climgi = climg EQ cl[i]
     limg   = TEMPORARY(limg) + LABEL_REGION(climgi, /ULONG,all_neighbors=all_beighbors) + climgi*mxlb
   ENDFOR
   progressbar->destroy

; re-label segments larger than min_size
   indices = where(histogram(limg,reverse_indices=r) ge min_size,count)
   limg=limg*0
   progressbar = Obj_New('progressbar', Color='blue', Text=' ',$
                 title='Re-labeling...',xsize=300,ysize=20)
   progressbar->start
   if count gt 0 then for i=1L,count-1 do begin
      if progressbar->CheckCancel() then begin
         progressbar->destroy
         return,-1
      endif
      pct = (i)*100.0/count
      progressbar->Update,fix(pct),text=strtrim(i)
      j = indices[i]       ;label of segment
      p = r[r[j]:r[j+1]-1] ;pixels in segment 
      limg[p]=i            ;new label
   endfor
   progressbar->destroy
   return, limg
end