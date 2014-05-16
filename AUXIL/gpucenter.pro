; docformat = 'rst'
; gpucenter.pro
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
;       center a kernel matrix::
;          Shawe-Taylor, J. and Cristianini, N. (2004).
;          Kernel Methods for Pattern Analysis.
;          Cambridge University Press. 
; :Params:
;       K_gpu:  in, required, type = GPU variable 
;          kernel matrix to be centered
; :Returns:
;       GPU variable, the centered kernel  
; :Requires: 
;       GPULIB                        
; :Author:
;       Mort Canty (2011)       
;-
function gpucenter, K_gpu
   m = (k_gpu.getdimensions())[0]
   ones_gpu = gpuMake_Array(m,m,value=1.0)
   tmp_gpu = ones_gpu##K_gpu + K_gpu##ones_gpu
   return,gpuAdd(1.0,K_gpu,-1.0/m,tmp_gpu,gpuTotal(K_gpu)/m^2)
end