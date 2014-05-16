; docformat = 'rst'
; gpuquadkernel_matrix.pro
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
;       Returns array of quadratic kernel functions 
;       for data matrices G1_gpu and G2_gpu. 
;       G1_gpu is an N x m array,
;       G2_gpu is an N x n array.
;       Returned array is n x m. 
;       
;       If G2_gpu not present then returns a symmetric, 
;       positive definite quadratic kernel matrix  
; :Returns: {GPUHANDLE}           
; :Params:
;       G1_gpu:  in, required, type={GPUHANDLE} 
;          data matrix
;       G2_gpu:  in, optional, type={GPUHANDLE}
;          data matrix       
; :Keywords:
;       bias: in, optional
;            bias (default 1)                    
; :Author:
;       Mort Canty (2009)
;       Juelich Research Center
;       m.canty@fz-juelich.de       
; :Uses:
;       GPULib
;-
function gpuquadkernel_matrix, G1_gpu, G2_gpu, BIAS=bias
   if n_params() eq 1 then G2_gpu=G1_gpu
   if n_elements(bias) eq 0 then bias = 1.0
   tmp_gpu = gpumatrix_multiply(G2_gpu,G1_gpu,/atranspose) 
   tmp_gpu = gpuadd(1.0,tmp_gpu,0.0,tmp_gpu,bias,lhs=tmp_gpu)
   K_gpu = gpumult(tmp_gpu,tmp_gpu)
   gpufree,tmp_gpu
   return,K_gpu
end   