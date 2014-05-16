; docformat = 'rst'
; gpugausskernel_matrix.pro
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
;       Returns array of Gaussian kernel functions 
;       for data matrices G1_gpu and G2_gpu. 
;       G1_gpu is an N x m array,
;       G2_gpu is an N x n array.       
;       If G2_gpu not present then returns a symmetric, 
;       positive definite Guassian kernel matrix  
; :Returns: 
;       {GPUHANDLE} to n x m array           
; :Params:
;       G1_gpu:  in, required, type={GPUHANDLE} 
;          data matrix
;       G2_gpu:  in, optional, type={GPUHANDLE}
;          data matrix       
; :Keywords:
;       gma: in,out,optional,type=float            
;            if not given, calculated from the data:
;            GMA=1/(2*(nscale*scale)^2) where
;            scale = average distance between observations 
;            in the feature space  
;       nscale: in, optional,type=float
;            multiple of scale for GMA (default 1.0)                             
; :Uses:
;       GPULib
; :Author:
;       Mort Canty (2009) 
;-
function gpugausskernel_matrix, G1_gpu, G2_gpu, $
                        GMA=gma, NSCALE=nscale
   COMPILE_OPT STRICTARR 
   if n_params() eq 1 then G2_gpu=G1_gpu
   if n_elements(nscale) eq 0 then nscale = 1.0
   m = (G1_gpu.getdimensions())[1]
   n = (G2_gpu.getdimensions())[1]  
   n1_gpu = gpumake_array(n,value=1.0)
   G12_gpu = gpumult(G1_gpu,G1_gpu)
   t1_gpu = gputotal(G12_gpu,1)
   gpuFree, G12_gpu
   K_gpu = gpumatrix_multiply(n1_gpu,t1_gpu,/btranspose)
   gpuFree, [t1_gpu,n1_gpu]
   G22_gpu = gpumult(G2_gpu,G2_gpu)
   t2_gpu = gputotal(G22_gpu,1)
   gpuFree, G22_gpu
   m1_gpu = gpumake_array(m,value=1.0)
   t3_gpu = gpumatrix_multiply(t2_gpu,m1_gpu,/btranspose)
   gpuFree, [t2_gpu,m1_gpu]
   K_gpu = gpuadd(K_gpu,t3_gpu,lhs=K_gpu)
   gpuFree, t3_gpu
   t4_gpu = gpumatrix_multiply(G2_gpu,G1_gpu,/atranspose)
   K_gpu = gpuadd(1.0,K_gpu,-2.0,t4_gpu,0.0,lhs=K_gpu)
   gpufree, t4_gpu
   if n_elements(gma) eq 0 then begin
      gpuAbs, K_gpu, Ka_gpu
      Ksqrt_gpu = gpusqrt(Ka_gpu)
      scale = gputotal(Ksqrt_gpu)/(m^2-m)
      gpufree,Ksqrt_gpu
      gpuFree,Ka_gpu
      gma = 1/(2*(nscale*scale)^2)
   endif 
   K_gpu = gpuexp(1.0,-gma,K_gpu,0.0,0.0,lhs=K_gpu)
   return, K_gpu   
end   