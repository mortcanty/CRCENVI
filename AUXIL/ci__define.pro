; docformat = 'rst'
; ci__define.pro
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
;       Object class to find thin closed contours in an
;       image band with combined Sobel-LoG filtering
; :Params:
;       image: in, required
;           image to be searched for contours
;       sigma: in, required
;           radius for LoG filter
; :Examples:         
;        ci = Obj_New("CI",image,sigma)                            
; :Author:
;       Mort Canty (2009)      
;-
Function CI::Init, image, sigma
; parameters --------------
   self.T1 = 32              ; low trace threshold
   self.T2 = 128             ; high trace threshold
   self.min_length = 50      ; minimum contour length
   self.max_length = 300     ; maximum contour length
   self.max_contours = 8000  ; maximum number of contours
; -------------------------
   sz = size(image)
   self.tile_count = 0
   self.num_cols = sz[1]
   self.num_rows = sz[2]
   self.num_tiles = sz[2]/50
   self.contour_image = Ptr_New(byte(image*0))
   image = float(image)
;Chain code difference lookup table
   for i=0,7 do self.cc_lut[i,*] = shift(indgen(8),i)
;First do a Sobel filter
   filt1 = fltarr(3,3)
   filt2 = fltarr(3,3)
   filt1[0,0]=1  & filt1[1,0]=2  & filt1[2,0]= 1
   filt1[0,1]=0  & filt1[1,1]=0  & filt1[2,1]= 0
   filt1[0,2]=-1 & filt1[1,2]=-2 & filt1[2,2]= -1
   filt2[0,0]=-1 & filt2[1,0]=0  & filt2[2,0]= 1
   filt2[0,1]=-2 & filt2[1,1]=0  & filt2[2,1]= 2
   filt2[0,2]=-1 & filt2[1,2]=0  & filt2[2,2]= 1
   edge_strength = bytscl(sqrt(convol(image,filt1,/center)^2+ $
                       convol(image,filt2,/center)^2))
;Then perform LoG filter
   r = round(3*sigma*sqrt(2))
   LoGfilt = fltarr(2*r+1,2*r+1)
   s2 = 2*sigma^2
   s6 = 2*!Pi*sigma^6
   for i=0,2*r do for j=0,2*r do begin
     xx = (i-r)^2+(j-r)^2
     LoGfilt[i,j] = (xx-s2)*exp(-xx/s2)/s6
   endfor
   temp = convol(image,LoGfilt,/center)
;Write edge strengths at zero-crossings to contour_image
   edges = where( (temp*shift(temp,1,0) lt 0) or (temp*shift(temp,0,1) lt 0) )
   (*self.contour_image)[edges] = edge_strength[edges]
;Normalize to +-4 sdev and store in contour_image
   mn = mean((*self.contour_image)[edges])
   sdev = sqrt(variance((*self.contour_image)[edges]))
   max = mn + 4*sdev
   min = mn - 4*sdev
   temp = ((*self.contour_image)[edges]-min)*255/(max-min)
   temp = temp<255
   (*self.contour_image)[edges] = byte(temp>0)
   Return, 1
End

Pro CI::Cleanup
   Ptr_Free, self.contour_image
End

Function CI::Get_Max_Contours
   return, self.max_contours
end

Function CI::Get_Max_Length
   return,self.max_length
end

Function CI::Get_Contour_Image
   return,*self.contour_image
End

Pro CI::Clear_Contour_Image
   *self.contour_image = (*self.contour_image)*0
End

Pro CI::Write_Contour, cont, grayscale
   pts = self->to_pixels(cont)
   for i=0,n_elements(pts[*,0])-1 do begin
       x= pts[i,0]
       y= pts[i,1]
      (*self.contour_image)[x,y]=grayscale
   endfor
End

Function CI::To_Pixels,cont
   result = intarr(cont.length,2)
   result[0,*] = cont.sp
   x = cont.sp[0]
   y = cont.sp[1]
   for i=1,cont.length-1 do begin
      case cont.code[i] of
         0: x=x+1
         1: begin
               x=x+1
               y=y-1
            end
         2: y=y-1
         3: begin
               x=x-1
               y=y-1
            end
         4: x=x-1
         5: begin
               x=x-1
               y=y+1
            end
         6: y=y+1
         7: begin
               x=x+1
               y=y+1
            end
      endcase
      result[i,*] = [x,y]
   endfor
   return,result
End

Function CI::To_Filtered_CC, cont
   kernel = [0.1,0.2,0.4,0.2,0.1]
   result = fltarr(n_elements(cont.code))
   result[0] = cont.code[0]
   for i=1, n_elements(cont.code)-1 do begin
      a = cont.code[i]
      b = a + 8
      if abs(a-result[i-1]) lt abs(b-result[i-1]) then $
         result[i]=a else result[i]=b
   endfor
   return, convol(result,kernel)
End

Function CI::To_Moments, cont
   hu = fltarr(8)
   pts = float(self->to_pixels(cont))
   m00 = n_elements(pts[*,0])
   m10 = total(pts[*,0])
   m01 = total(pts[*,1])
   xbar = m10/m00
   ybar = m01/m00
   m11 = total(pts[*,0]*pts[*,1])
   m20 = total(pts[*,0]^2)
   m02 = total(pts[*,1]^2)
   m21 = total(pts[*,0]*pts[*,0]*pts[*,1])
   m12 = total(pts[*,0]*pts[*,1]*pts[*,1])
   m30 = total(pts[*,0]^3)
   m03 = total(pts[*,1]^3)
   mu00 = m00
   eta11 = (m11 - xbar*m01)/mu00^2
   eta20 = (m20 - xbar*m10)/mu00^2
   eta02 = (m02 - ybar*m01)/mu00^2
   eta21 = (m21 - 2*xbar*m11 - ybar*m20 + 2*xbar*xbar*m01)/mu00^(5.0/2.0)
   eta12 = (m12 - 2*ybar*m11 - xbar*m02 + 2*ybar*ybar*m10)/mu00^(5.0/2.0)
   eta30 = (m30-3*xbar*m20 + 2*xbar*xbar*m10)/mu00^(5.0/2.0)
   eta03 = (m03-3*ybar*m02 + 2*ybar*ybar*m01)/mu00^(5.0/2.0)
   hu[0] = mu00
   hu[1] = eta20 + eta02
   hu[2] = (eta20 - eta02)^2 + 4*eta11^2
   hu[3] = (eta30 - 3*eta12)^2 + (3*eta21 - eta03)^2
   hu[4] = (eta30 + eta12)^2 + (eta21 + eta03)^2
   return, alog(hu)
End


Function CI::Trace_Contour
;Basic contour structure
   c = { sp:      intarr(2),                 $   ; starting point coordinates
         length:  0L,                        $   ; number of pixels
         closed:  -1L,                       $   ; set to zero while tracing, one when closed, -1 when no more contours
         code:    bytarr(self.max_length),   $   ; chain code
         icode:   bytarr(self.max_length)    }   ; rotationally invariant (difference) chain code
   pts = intarr(self.max_length+1,2)
;Get a starting point
   addr = where((*self.contour_image)[*,self.tile_count*self.num_rows/self.num_tiles: $
                                       (self.tile_count+1)*self.num_rows/self.num_tiles-1] gt self.T2, count)
   if (count le 1) then begin
;No more starting points in this tile
;Bump the tile counter if not processing last tile and return with c.closed=0
;Otherwise return with c.closed=-1 to indicate no more starting points
      if self.tile_count lt self.num_tiles-1 then begin
         c.closed=0
         self.tile_count=self.tile_count+1
      endif
      return, c
   endif
;Get the xy address of the starting point
   x = addr[0] mod self.num_cols
   y = addr[0]/self.num_cols + self.tile_count*self.num_rows/self.num_tiles
   c.sp = [x,y]
   pts[0,*]=[x,y]
;Delete it
   (*self.contour_image)[x,y] = 0
;Trace contour, deleting each pixel after it is found
   c.closed = 0
   done = 0
   while (done eq 0) and (c.closed eq 0) do begin
      index = where((*self.contour_image)[x-1:x+1,y-1:y+1] gt self.T1, count)
;    no more valid pixels, so done
      if count eq 0 then done = 1 $
;    otherwise continue tracing
      else begin
         case index[0] of
         0: begin
               c.code[c.length] = 3 ; northwest
               x = x-1
               y = y-1
            end
         1: begin
               c.code[c.length] = 2 ; north
               y = y-1
            end
         2: begin
               c.code[c.length] = 1 ; northeast
               x = x+1
               y = y-1
            end
         3: begin
               c.code[c.length] = 4 ; west
               x = x-1
            end
         5: begin
               c.code[c.length] = 0 ; east
               x = x+1
            end
         6: begin
               c.code[c.length] = 5 ; southwest
               x = x-1
               y = y+1
            end
         7: begin
               c.code[c.length] = 6 ; south
               y = y+1
            end
         8: begin
               c.code[c.length] = 7 ; southeast
               x = x+1
               y = y+1
            end
         endcase
;       delete pixel
         (*self.contour_image)[x,y] = 0
         c.length=c.length+1
         pts[c.length,*]=[x,y]
      endelse
;Check that contour is not too large
      if (c.length eq self.max_length-1) then done=1
;Check for (almost) closed contour
      if (c.length gt self.min_length) and (done eq 0) then begin
         temp = [[pts[0:c.length-self.min_length,0]-x],[pts[0:c.length-self.min_length,1]-y]]
         index = where(sqrt(total(temp^2,2)) le 3, count)
         if count gt 0 then begin
            c.closed = 1
            c.sp = pts[index[0],*]
            c.code[0:c.length-index[0]] = c.code[index[0]:c.length]
            c.length = c.length-index[0]+2  ; contour length is chain code length + 1
;Add the rotationally invariant chain code
            for i=0,c.length-1 do c.icode[i] = self.cc_lut[c.code[i],c.code[(i+1) mod c.length]]
         endif
      endif
   endwhile
   return, c
End

Pro CI__Define
class   = { CI, $
           contour_image: Ptr_New(), $
           cc_lut: intarr(8,8),$
           min_length:    0L,  $
           max_length:    0L,  $
           T1:            0L,  $
           T2:            0L,  $
           max_contours:  0L,  $
           tile_count:    0L,  $
           num_tiles:     0L,  $
           num_rows:      0L,  $
           num_cols:      0L }
End