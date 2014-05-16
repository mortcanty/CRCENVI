; docformat = 'rst'
; hu_moments.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

function SM,A,p,q
; calculate discrete moments of A
   sz = size(A)
   x = findgen(sz[1],sz[2]) mod sz[1]
   y = findgen(sz[2])
   return, y^q##transpose(total(x^p*A,1,/double))
end

;+
; :Description:
;    Takes a rectangular array A of gray values and returns
;    the 7 Hu invariant moments      
; :Params:
;      A:  in, required 
;         array of gray-scale values            
; : Keywords:
;     log: in, optional
;         natural logarithms of the moments
;         are returned if set                        
; :Author:
;      Mort Canty (2009)      
;-
function hu_moments, A, log=log
   if n_elements(log) eq 0 then log=0
; standard moments, corrected after Liao and Pawlak
   dx = 1./((size(A))[1]-1)
   dy = 1./((size(A))[2]-1)
   M00=SM(A,0,0)
   M01=SM(A,0,1)
   M10=SM(A,1,0)
   M20=SM(A,2,0)+dx*dx*M00/12
   M02=SM(A,0,2)+dy*dy*M00/12
   M11=SM(A,1,1)
   M30=SM(A,3,0)+dx*dx*M10/4
   M21=SM(A,2,1)+dx*dx*M01/12
   M12=SM(A,1,2)+dy*dy*M10/12
   M03=SM(A,0,3)+dy*dy*M01/4
; normalized central moments
   Xbar=M10/M00
   Ybar=M01/M00
   Eta11=(M11-Ybar*M10)/M00^2
   Eta20=(M20-Xbar*M10)/M00^2
   Eta02=(M02-Ybar*M01)/M00^2
   Eta30=(M30-3*Xbar*M20+2*Xbar*Xbar*M10)/M00^2.5
   Eta03=(M03-3*Ybar*M02+2*Ybar*Ybar*M01)/M00^2.5
   Eta21=(M21-2*Xbar*M11-Ybar*M20+2*Xbar*Xbar*M01)/M00^2.5
   Eta12=(M12-2*Ybar*M11-Xbar*M02+2*Ybar*Ybar*M10)/M00^2.5
; Hu invariant moments
   moms = [Eta20+Eta02, $
           (Eta20-Eta02)^2+4*Eta11^2, $
           (Eta30-3*Eta12)^2+(3*Eta21-Eta03)^2, $
           (Eta30+Eta12)^2+(Eta21+Eta03)^2, $
           (Eta30-3*Eta12)*(Eta30+Eta12)*((Eta30+Eta12)^2-3*(Eta21+Eta03)^2)+(3*Eta21-Eta03)*(Eta21+Eta03)*(3*(Eta30+Eta12)^2-(Eta21+Eta03)^2), $
           (Eta20-Eta02)*((Eta30+Eta12)^2-(Eta21+Eta03)^2)+4*Eta11*(Eta30+Eta12)*(Eta21+Eta03), $
           (3*Eta21-Eta03)*(Eta30+Eta12)*((Eta30+Eta12)^2-3*(Eta21+Eta03)^2)+(3*Eta12-Eta30)*(Eta21+Eta03)*(3*(Eta30+Eta12)^2-(Eta21+Eta03)^2)]
   if log then return, alog(abs(moms)+1e-30) else return, moms                  
end