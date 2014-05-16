; docformat = 'rst'
; polykernel_matrix.pro
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
;       Returns array of polynomial kernel functions 
;       for data matrices G1 and G2. 
;       G1 is an N x m array,
;       G2 is an N x n array.
;       Returned array is n x m. 
;       If G2 not present then returns a symmetric, 
;       positive definite kernel matrix       
; :Params:
;       G1:  in, required 
;          data matrix
;       G2:  in, optional  
;          data matrix       
; :Keywords:
;       degree: in, optional
;            degree of polynomial (default 2)
;       bias: in, optional
;            bias (default 1)                       
; :Author:
;       Mort Canty (2009)
;       Juelich Research Center
;       m.canty@fz-juelich.de       
; :Uses:
;       GPULib
;-
function polykernel_matrix, G1, G2, DEGREE=degree, BIAS=bias
   if n_params() eq 1 then G2 = G1
   if n_elements(degree) eq 0 then degree = 2.0
   if n_elements(bias) eq 0 then bias = 1.0
   return, (G1##transpose(G2)+bias)^degree 
end   