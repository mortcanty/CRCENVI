; docformat = 'rst'
; gpukernel_matrix.pro
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
;       Returns array of kernel functions 
;       for data matrices G1 and G2. 
;       G1 is an N x m array,
;       G2 is an N x n array.
;       Returned array is n x m. 
;       If G2 not present then returns a symmetric, 
;       positive definite kernel matrix       
; :Params:
;       G1_gpu:  in, required, type={GPUHANDLE} 
;          data matrix
;       G2_gpu:  in, optional, type={GPUHANDLE}
;          data matrix     
; :Keywords:
;       gma: in,out,optional,type=float            
;            if not given and KERNEL=1, calculated from the data:
;            GMA=1/(2*(nscale*scale)^2) where scale = average 
;            distance between observations in the input space,
;            otherwise defaults to 1/N 
;       nscale: in, optional,type=float
;            multiple of scale for GMA when KERNEL=1 (default 1.0)    
;       kernel: in, optional, type=integer
;            the kernel used
;            0: linear
;            1: Gaussian (default)
;            2: polynomial
;            3: sigmoid
;       degree: in, optional, type=integer
;            degree of pylynomial kernel (default 2)
;       bias: in, optional, type=float 
;            bias of polynomial or sigmoid kernel (default 1.0)      
; :Requires: 
;       GPULIB                                                      
; :Author:
;       Mort Canty (2009)      
;-
function gpukernel_matrix,G1_gpu,G2_gpu,KERNEL=kernel,GMA=gma,NSCALE=nscale,DEGREE=degree,BIAS=bias
COMPILE_OPT STRICTARR
   if n_params() eq 1 then G2_gpu = G1_gpu
   if n_elements(nscale) eq 0 then nscale = 1.0
   if n_elements(kernel) eq 0 then kernel = 1
   if n_elements(degree) eq 0 then degree = 2
   if n_elements(bias) eq 0 then bias = 1.0
   p = (G1_gpu.getdimensions())[0]
   m = (G1_gpu.getdimensions())[1]
   n = (G2_gpu.getdimensions())[1]  
   tp = G1_gpu.getType()  
   case kernel of
      0: begin
            K_gpu = gpumatrix_multiply(G2_gpu,G1_gpu,/atranspose)
            gma = 0.0
         end   
      1: begin  
            dbl = (tp eq 5) ? 1 : 0         
            n1_gpu = gpumake_array(n,value=1.0,double=dbl)                    
            K_gpu = gpumatrix_multiply(n1_gpu,gpuTotal(G1_gpu*G1_gpu,1),/btranspose)                                        
            m1_gpu = gpumake_array(1,m,value=1.0,double=dbl)
            K_gpu = K_gpu + m1_gpu##gputotal(G2_gpu*G2_gpu,1)                           
            K_gpu = gpuadd(1.0,K_gpu,-2.0,G1_gpu##gpuTranspose(G2_gpu),0.0)
            if n_elements(gma) eq 0 then begin                     
               scale = gputotal(gpuSqrt(gpuAbs(K_gpu)))/(m^2-m)
               gma = 1/(2*(nscale*scale)^2)
            endif   
            K_gpu = gpuexp(1.0,-gma,K_gpu,0.0,0.0)  
         end
      2: begin
            if n_elements(gma) eq 0 then gma=1.0/p        
            K1_gpu = g1_gpu##gpuTranspose(G2_gpu)           
            K1_gpu = gpuadd(gma,K1_gpu,0.0,K1_gpu,bias,lhs=K1_gpu)
            K_gpu = gpuCopy(K1_gpu)
            for i=1,degree-1 do K_gpu = K_gpu*K1_gpu
         end
      3: begin            
            K_gpu = G1_gpu##gpuTranspose(G2_gpu)
            if n_elements(gma) eq 0 then $
               gma = 1/mean(abs(gpuGetArr(K_gpu)))
            K_gpu = gpuadd(gma,K_gpu,0.0,K_gpu,bias)
            K_gpu = gpusinh(K_gpu)/gpucosh(K_gpu) 
         end
   endcase    
   return, K_gpu          
end   