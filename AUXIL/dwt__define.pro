; docformat = 'rst'
; dwt__define.pro
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
;       Discrete wavelet transform object class using
;       Daubechies wavelets for construction of pyramid
;       representations of images::
;         Ranchin, T. and Wald, L. (2000). 
;         Fusion of high spatial and spectral resolution
;         images: the ARSIS concept and its implementation. 
;         Photogrammetric Engineering and Remote Sensing, 
;         66(1), 49–61.
; :Params:
;       image: in,required
;          grayscale image to be transformed
; :Examples:         
;        dwt = Obj_New("DWT",image)  
; :Uses:
;       PHASE_CORR::
;       WARP_SHIFT                                
; :Author:
;       Mort Canty (2009)      
;-
Function DWT::Init, image
; D4 refinement coefficients
   self.D4[0] = (1+sqrt(3))/8
   self.D4[1]= (3+sqrt(3))/8
   self.D4[2]= (3-sqrt(3))/8
   self.D4[3] = (1-sqrt(3))/8
; D6 refinement coefficients
   self.D6[0] = (1+sqrt(10)+sqrt(5+2*sqrt(10)))/32
   self.D6[1] = (5+sqrt(10)+3*sqrt(5+2*sqrt(10)))/32
   self.D6[2] = (10-2*sqrt(10)+2*sqrt(5+2*sqrt(10)))/32
   self.D6[3] = (10-2*sqrt(10)-2*sqrt(5+2*sqrt(10)))/32
   self.D6[4] = (5+sqrt(10)-3*sqrt(5+2*sqrt(10)))/32
   self.D6[5] = (1+sqrt(10)-sqrt(5+2*sqrt(10)))/32
; D8 refinement coefficients
   self.D8[0] =  0.2304/sqrt(2)
   self.D8[1] =  0.7148/sqrt(2)
   self.D8[2] =  0.6309/sqrt(2)
   self.D8[3] = -0.0280/sqrt(2)
   self.D8[4] = -0.1870/sqrt(2)
   self.D8[5] =  0.0308/sqrt(2)
   self.D8[6] =  0.0329/sqrt(2)
   self.D8[7] = -0.0106/sqrt(2)
; default is D4
   self->set_coeff, 4
   self.compressions = 0
   self.max_compressions = 4
   self.sz = size(image)
; ignore image edges if dimension not
; divisible by 2^max_compressions
   r = 2^self.max_compressions
   self.num_cols = r*(self.sz[1]/r)
   self.num_rows = r*(self.sz[2]/r)
   self.image = Ptr_New(float(image[0:self.num_cols-1,0:self.num_rows-1]))
   Return, 1
End

Pro DWT::Cleanup
   Ptr_Free, self.image
   Ptr_Free, self.H
   Ptr_Free, self.G
End

Pro DWT::Set_Coeff, D
   Ptr_Free, self.H
   Ptr_Free, self.G
   case D of
      6: begin
            self.H = Ptr_New(reverse(self.D6))
            self.G = Ptr_New(self.D6*[-1,1,-1,1,-1,1])
            self.D = 6
         end
      8: begin
            self.H = Ptr_New(reverse(self.D8))
            self.G = Ptr_New(self.D8*[-1,1,-1,1,-1,1,-1,1])
            self.D = 8
         end
     else: begin ; defaults to D4
            self.H = Ptr_New(reverse(self.D4))
            self.G = Ptr_New(self.D4*[-1,1,-1,1])
            self.D = 4
           end
   endcase
end

Pro DWT::Show_Image,wn
   order = !order
   !order = 1
   window,wn,xsize=self.num_cols,ysize=self.num_rows
   wset,wn
   tv, bytscl(*self.image)
   !order = order
End

Pro DWT::Inject, im, quadrant=quadrant
if n_elements(quadrant) eq 0 then quadrant = 0
; overwrite quadrant
   n = self.num_cols/2^self.compressions
   m = self.num_rows/2^self.compressions
   case quadrant of
      0: (*self.image)[0:n-1,0:m-1]     = im[0:n-1,0:m-1]
      1: (*self.image)[n:2*n-1,0:m-1]   = im[0:n-1,0:m-1]
      2: (*self.image)[0:n-1,m:2*m-1]   = im[0:n-1,0:m-1]
      3: (*self.image)[n:2*n-1,m:2*m-1] = im[0:n-1,0:m-1]
   endcase
End

Pro DWT::InjectExact, im, offset=offset
   n = self.num_cols/2^self.compressions
   m = self.num_rows/2^self.compressions
   old = (*self.image)[0:n-1,0:m-1]
   new = im[0:n-1,0:m-1]
   offset = phase_corr(new,old,/subpixel,display=10)
   (*self.image)[0:n-1,0:m-1] = warp_shift(new,offset)
End

Pro DWT::Set_Compressions,cm
   self.compressions = cm
End

Function DWT::Get_Compressions
   return, self.compressions
End

Function DWT::Get_Num_Cols
   return, self.num_cols/2^(self.compressions)
End

Function DWT::Get_Num_Rows
   return, self.num_rows/2^(self.compressions)
End

Function DWT::Get_Image
   im = fltarr(self.sz[1],self.sz[2])
   im[0:self.num_cols-1,0:self.num_rows-1]=*self.image
   return, im
End

Function DWT::Get_Quadrant, k, no_edges=no_edges
   if n_elements(no_edges) eq 0 then no_edges=0
; get compressed image (as 2D array) or innermost wavelet coefficients as vector
   n = self.num_cols/2^(self.compressions)
   m = self.num_rows/2^(self.compressions)
   if k eq 0 then f = (*self.image)[0:n-1,0:m-1] $
   else begin
      f = fltarr(m*n)
      case k of
        1: f = (*self.image)[n:2*n-1,0:m-1]
        2: f = (*self.image)[0:n-1,m:2*m-1]
        3: f = (*self.image)[n:2*n-1,m:2*m-1]
      endcase
   endelse
   if no_edges then f = f[2:n-3,2:m-3]
   return, f[*]
End

Pro DWT::Normalize_WC, a, b
; normalize wavelet coefficients at all levels
   if self.compressions gt 0 then for c=1,self.compressions do begin
      n = self.num_cols/2^c
      m = self.num_rows/2^c
      fg = a[0]*(*self.image)[n:2*n-1,0:m-1]+b[0]
      (*self.image)[n:2*n-1,0:m-1] = fg
      gf= a[1]*(*self.image)[0:n-1,m:2*m-1]+b[1]
      (*self.image)[0:n-1,m:2*m-1] = gf
      gg= a[2]*(*self.image)[n:2*n-1,m:2*m-1]+b[2]
      (*self.image)[n:2*n-1,m:2*m-1] = gg
   endfor
End

Pro DWT::Filter
if self.compressions lt self.max_compressions $
then begin
; single application of filter bank
   H = *self.H & G = *self.G
   m = self.num_rows/2^(self.compressions)
   n = self.num_cols/2^(self.compressions)
   f0 = (*self.image)[0:n-1,0:m-1]  
; temporary arrays for wavelet coefficients
   f1 = fltarr(n,m/2)
   g1 = f1*0
   ff1 = fltarr(n/2,m/2)
   fg1 = ff1*0 & gf1 = ff1*0 & gg1 = ff1*0  
; filter columns and downsample
   ds = indgen(m/2)*2+1
   for i=0,n-1 do begin
      temp = convol(transpose(f0[i,*]),H,/edge_wrap)
      f1[i,*] = temp[ds]
      temp = convol(transpose(f0[i,*]),G,/edge_wrap)
      g1[i,*] = temp[ds]
   endfor   
; filter rows and downsample
   ds = indgen(n/2)*2+1
   for i=0,m/2-1 do begin
      temp = convol(f1[*,i],H,/edge_wrap)
      ff1[*,i] = temp[ds]
      temp = convol(f1[*,i],G,/edge_wrap)
      fg1[*,i] = temp[ds]
      temp = convol(g1[*,i],H,/edge_wrap)
      gf1[*,i] = temp[ds]
      temp = convol(g1[*,i],G,/edge_wrap)
      gg1[*,i] = temp[ds]
   endfor
   f0[0:n/2-1,0:m/2-1] = ff1[*,*]
   f0[n/2:n-1,0:m/2-1] = fg1[*,*]
   f0[0:n/2-1,m/2:m-1] = gf1[*,*]
   f0[n/2:n-1,m/2:m-1] = gg1[*,*]
   (*self.image)[0:n-1,0:m-1] = f0
   self.compressions = self.compressions+1
endif
End

Pro DWT::Invert
if self.compressions gt 0 then begin
; single application of inverse filter bank
   H = *self.H
   G = *self.G
   H = reverse(H)
   G = reverse(G)
   m = self.num_rows/2^(self.compressions-1)
   n = self.num_cols/2^(self.compressions-1)
; temporary arrays for wavelet coefficients
   f0  = (*self.image)[0:n-1,0:m-1]
   ff1 = f0[0:n/2-1,0:m/2-1]
   fg1 = f0[n/2:n-1,0:m/2-1]
   gf1 = f0[0:n/2-1,m/2:m-1]
   gg1 = f0[n/2:n-1,m/2:m-1]
   f1 = fltarr(n,m/2)
   g1 = fltarr(n,m/2)
; upsample and filter the rows
   for i=0,m/2-1 do begin
      a = reform(transpose([[ff1[*,i]],[fltarr(n/2)]]),n)
      b = reform(transpose([[fg1[*,i]],[fltarr(n/2)]]),n)
      f1[*,i] = convol(a,H,/edge_wrap) + convol(b,G,/edge_wrap)
   endfor
   for i=0,m/2-1 do begin
      a = reform(transpose([[gf1[*,i]],[fltarr(n/2)]]),n)
      b = reform(transpose([[gg1[*,i]],[fltarr(n/2)]]),n)
      g1[*,i] = convol(a,H,/edge_wrap) + convol(b,G,/edge_wrap)
   endfor
; upsample and filter the columns
   for i=0,n-1 do begin
      a = reform(transpose([[transpose(f1[i,*])],[fltarr(m/2)]]),m)
      b = reform(transpose([[transpose(g1[i,*])],[fltarr(m/2)]]),m)
      f0[i,*] = 4*( convol(a,H,/edge_wrap) + convol(b,G,/edge_wrap) )
   endfor
   (*self.image)[0:n-1,0:m-1] = f0
   self.compressions = self.compressions-1
endif
End

Pro DWT__Define
class =  { DWT, $
           image: Ptr_New(), $
           H: Ptr_New(),     $
           G: Ptr_New(),     $
           compressions: 0,  $
           max_compressions: 0,$
           D: 0,          $
           sz: intarr(5), $
           D4: fltarr(4), $
           D6: fltarr(6), $
           D8: fltarr(8), $
           num_rows: 0L,  $
           num_cols: 0L }
End

