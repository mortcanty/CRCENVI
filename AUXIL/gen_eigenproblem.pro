; docformat = 'rst'
; gen_eigenproblem.pro 
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
;       Solve the generalized eigenproblem 
;       C##a = lambda*B##a using Cholesky factorization   
; :Params:
;       C:  in, required 
;          symmetric matrix
;       B:  in, required   
;          symmetric positive definite matrix 
;       A:  out, required
;          matrix of eigenvectors (columns)
;       lambda: out, required
;           eigenvalues, largest to smallest                         
; :Author:
;       Mort Canty (2009)    
;-
pro gen_eigenproblem, C,B,A,lambda
; solve the generalized eigenproblem C##a = lambda*B##a
   choldc, B, P, /double
   for i=1L,(size(B))[1]-1 do B[i,0:i]=[fltarr(i),P[i]]
   B[0,0]=P[0]
   Li  = invert(B,/double)
   D = Li ## C ## transpose(Li)
; ensure symmetry after roundoff errors
   D = (D+transpose(D))/2
   lambda = eigenql(D,/double,eigenvectors=A)
; eigenvectors are in columns of A
   A = transpose(A##Li)
end