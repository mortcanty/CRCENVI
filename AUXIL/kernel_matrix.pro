; docformat = 'rst'
; kernel_matrix.pro
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
;       G1:  in, required 
;          data matrix
;       G2:  in, optional  
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
; :Author:
;       Mort Canty (2009)      
;-
function kernel_matrix,G1,G2,KERNEL=kernel,GMA=gma,NSCALE=nscale,DEGREE=degree,BIAS=bias
   if n_params() eq 1 then G2 = G1
   if n_elements(nscale) eq 0 then nscale = 1.0
   if n_elements(kernel) eq 0 then kernel = 1
   if n_elements(degree) eq 0 then degree = 2
   if n_elements(bias) eq 0 then bias = 1.0
   case kernel of
      0: begin
            gma = 0.0
            return, G1##transpose(G2)
         end   
      1: begin   
            m = n_elements(G1[0,*]) 
            n = n_elements(G2[0,*])
            K = transpose(total(G1^2,1))##(intarr(n)+1)
            K = K + transpose(intarr(m)+1)##total(G2^2,1)
            K = K - 2*G1##transpose(G2)            
            if n_elements(gma) eq 0 then begin
               scale = total(sqrt(abs(K)))/(m^2-m)
               gma = 1/(2*(nscale*scale)^2)
            endif  
            return, exp(-gma*K)
         end
       2: begin
             if n_elements(gma) eq 0 then gma=1.0/n_elements(G1[*,0])
             return, (gma*G1##transpose(G2)+bias)^degree 
          end
       3: begin
             K = G1##transpose(G2)
             if n_elements(gma) eq 0 then gma = 1/mean(abs(K))
             return, tanh(gma*K+bias)
          end
   endcase               
end   