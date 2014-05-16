; docformat = 'rst'
; phase_corr.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details. 
;+
; :Description:
;       Returns relative offset [xoff,yoff] of two images
;       using phase correlation.
;       Maximum offset should not exceed +- 5 pixels
;       in each dimension
;       Returns -1 if dimensions are not equal
;       
;       Ref: H, Shekarforoush et al. (1995) INRIA 2707
; :Params:
;       im1: input,required 
;       im2: input,required 
; :Keywords:      
;       display:  input,optional
;          show a surface plot if the
;          correlation in window with display
;          number display
;       subpixel: input, optional
;          returns result to subpixel accuracy if
;          set, otherwise nearest integer (default)  
; :Author:
;       Mort Canty (2009) 
;-
FUNCTION phase_corr, im1, im2, display = display, subpixel = subpixel

if n_elements(subpixel) eq 0 then subpixel = 0
if n_elements(display) eq 0 then display = 0
sz1 = size(im1)
sz2 = size(im2)
if (sz1[1] eq sz2[1]) and (sz1[2] eq sz2[2]) then begin
   f1 = fft(im1, /double)
   f2 = fft(im2, /double)
   den = abs(f2*conj(f1))
   indices = where(den eq 0, count)
   if count gt 0 then den[indices]= 0.0000001
; back-transformed cross power spectrum
   g =  abs(fft( f2*conj(f1)/den, /inverse, /double ))
; shift 10 pixels
   g = shift(g,[10,10])
; get maximum position
   _ = max(g,pos)  
   xoff = pos mod sz1[1]
   yoff = pos/sz1[1]
; get 5x5 window centered at maximum
   gg = g[xoff-2:xoff+2,yoff-2:yoff+2]
; calculate offsets to nearest integer
   xoff=xoff-10.0
   yoff=yoff-10.0
   
   if subpixel then begin
; subtract noise level from window
      nlev = total(g[20:sz1[1]-1,20:sz1[2]-1])/((sz1[1]-20)*(sz1[2]-20))
      gg = gg - nlev
; get center of gravity in window
      arr = fltarr(5,5)
      for i=0,4 do arr[*,i] = indgen(5)
; correct offset to subpixel accuracy
      xoff = xoff+total(gg*arr)/total(gg)-2.0
      yoff = yoff+total(gg*transpose(arr))/total(gg)-2.0
   end
   
   if display then begin
      isurface,g[0:20,0:20],view_title='Phase correlation'
      if subpixel then begin
         isurface,gg
      endif   
   endif
   
   return, [xoff,yoff]
end else return, -1
END