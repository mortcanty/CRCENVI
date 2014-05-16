 ; docformat = 'rst'
; atwt__define.pro
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
;       A trous wavelet transform object class using
;       cubic spline wavelet for panchromatic
;       sharpening of multispectral images::
;           Aiazzi, B., Alparone, L., Baronti, S., and 
;           Garzelli, A. (2002). Context-driven
;           fusion of high spatial and spectral resolution 
;           images based on oversampled ;multiresolution 
;           analysis. IEEE Transactions on Geoscience and 
;           Remote ;Sensing, 40(10), 2300–2312.
; :Params:
;       image: in,required
;          grayscale image to be transformed
; :Examples:         
;        dwt = Obj_New("ATWT",image)                            
; :Author:
;       Mort Canty (2009)      
;-
Function ATWT::Init, image
   self.transforms = 0
; cubic spline filter
   self.H[0] = 1.0/16
   self.H[1] = 1.0/4
   self.H[2] = 3.0/8
   self.H[3] = 1.0/4
   self.H[4] = 1.0/16
   sz = size(image)
   self.num_cols = sz[1]
   self.num_rows = sz[2]
   self.images = PtrArr(5)
   self.images[0] = Ptr_New(float(image))
   Return, 1
End

Pro ATWT::Cleanup
   Ptr_Free, self.images
End

Pro ATWT::Show_Image,i,wn
   order = !order
   !order = 1
   window,wn,xsize=self.num_cols,ysize=self.num_rows
   wset,wn
   tv, bytscl(*self.images[i])
   !order = order
End

Pro ATWT::Inject, im
; overwrite lowest level
   n = self.num_cols
   m = self.num_rows
   *self.images[0] = im[0:n-1,0:m-1]
End

Pro ATWT::Set_Transforms, ts
   self.transforms = ts
End

Function ATWT::Get_Transforms
   return, self.transforms
End

Function ATWT::Get_Num_Cols
   return, self.num_cols
End

Function ATWT::Get_Num_Rows
   return, self.num_rows
End

Function ATWT::Get_Image,i
   return, *self.images[i]
End

Pro ATWT::Normalize_WC, a, b
; normalize wavelet coefficients at all levels
if self.transforms gt 0 then for i=1,self.transforms do $
      *self.images[i] = a*(*self.images[i])+b
End

Pro ATWT::Compress
if self.transforms lt 4 then begin
; upsample the filter
   self.transforms = self.transforms + 1
   n = 5*(2^(self.transforms-1)+1)
   H = reform(transpose([[self.H],[fltarr(5,2^(self.transforms-1))]]),n)
   n1 = n_elements(H)-2^(self.transforms-1)
   H = H[0:n1-1]
 ; temporary array
   f1 = fltarr(self.num_cols,self.num_rows)
   ff1 = f1*0
; filter columns
   f0 = *self.images[0]
   for i=0,self.num_cols-1 do $
      f1[i,*] = convol(transpose(f0[i,*]),H,/edge_wrap)
; filter rows
   for j=0,self.num_rows-1 do $
      ff1[*,j] = convol(f1[*,j],H,/edge_wrap)
   self.images[self.transforms] = Ptr_New(*self.images[0]-ff1)
   *self.images[0] = ff1
endif
End

Pro ATWT::Expand
if self.transforms gt 0 then begin
   *self.images[0] = *self.images[0]+*self.images[self.transforms]
   Ptr_Free, self.images[self.transforms]
   self.transforms = self.transforms - 1
end
End

Pro ATWT__Define
class =  { ATWT, $
           images: PtrArr(5), $
           H: fltarr(5), $
           transforms: 0, $
           sz: intarr(5),$
           num_rows: 0L, $
           num_cols: 0L }
End

