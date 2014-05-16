; docformat = 'rst'
; rule_convert.pro
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
;       Convert an ENVI rule image to class
;       membership probabilities
;  :Uses:
;       ENVI     
; :Author:
;      Mort Canty (2009) 
;-
pro rule_convert

envi_select, title='Choose rule image', $
    fid=fid, dims=dims,pos=pos,/no_dims,/no_spec
if (fid eq -1) then begin
   print, 'cancelled'
   return
endif

map_info = envi_get_map_info(fid=fid)

num_cols = dims[2]-dims[1]+1
num_rows = dims[4]-dims[3]+1
num_classes = n_elements(pos)
num_pixels = num_cols*num_rows

; get MaxLike rule image
rule_image = dblarr(num_pixels, num_classes)
for i=0,num_classes-1 do rule_image[*,i] = $
    envi_get_data(fid=fid,dims=dims,pos=pos[i])

; exponentiate and normalize
prob_image = exp(rule_image - alog(1.0/num_classes))
den = total(prob_image,2,/double)
for i=0,num_classes-1 do  prob_image[*,i]= $
    prob_image[*,i]/den

; write to memory
envi_enter_data, $
  reform(bytscl(prob_image,min=0.0,max=1.0),num_cols, $
   num_rows,num_classes),map_info=map_info

end

