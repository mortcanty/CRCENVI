; docformat = 'rst'
; mean_shift.pro
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
;       Determines the mode of the indexed vector.
;       The vectors are in the common block 
;       variable DATA. The bandwidth is in the
;       common block variable HS
; :Params:
;        idx: in, required  
;           pixel index
;        cpts: out, required  
;           labeled pixels within distance HS/2 
;           along path to mode and
;           within distance HS of mode
;        cpts_max: out, required 
;           maximum labeled pixel                    
; :Author:
;      Mort Canty (2009) 
;-
function mean_shift,idx,cpts,cpts_max
   common cblock, data,hs,nc,nr,nb,n
   cpts=bytarr(nc,nr)
; initialize mean
   i = idx mod nc
   j = idx/nc
   m  = (data[i,j,*])[*]
   dm = 100.0
   iter=0L
   cpts_max=0L
; iterate
   while dm gt 0 and iter lt 100 do begin
      bi = i-hs>0
      ei = i+hs<nc-1
      bj = j-hs>0
      ej = j+hs<nr-1
      dta = data[bi:ei,bj:ej,*]
      nd = n_elements(dta)/(nb+2)
      dta = reform(dta,nd,nb+2)
      ones = fltarr(nd)+1
      m1 = m
      ms = transpose(ones##m)
      d2 = total((dta-ms)^2,2)
      indices = where( d2 le hs^2, count )
      if count gt 0 then begin
         ii = indices mod (ei-bi+1)
         jj = indices/(ei-bi+1)
         cpts_max = cpts_max > ( ((bj+jj)*nc+bi+ii)[count-1]+1  < n-1 )
;       update mean
         m = round( total(dta[indices,*],1)/count )
      endif
;   flag pixels near the current path      
      indices = where( d2 le hs^2/4, count )
      if count gt 0 then begin
         ii = indices mod (ei-bi+1)
         jj = indices/(ei-bi+1)
         cpts[bi+ii,bj+jj]=1B 
      endif          
      i = m[nb]
      j = m[nb+1]
      dm = max(abs(m-m1))
      iter++
   endwhile
   return, m
end