
; docformat = 'rst'
; gpuConvol2D.pro
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
;       Convolution of a two dimensional array
;       on the GPU in the frequency domain                
; :Params:
;       a_gpu:  in, required, type={GPUHANDLE} 
;          array to be convolved, overwritten with convolution
;       k_gpu:  in, required, type={GPUHANDLE}
;          convolution kernel (padded and centered)      
; :Keywords:
;                          
; :Uses:
;       GPULib
; :Author:
;       Mort Canty (20013) 
;-
function gpuConvol2D, a_gpu, k_gpu

   COMPILE_OPT IDL2

   a_cols = (a_gpu.getdimensions())[0]
   a_rows = (a_gpu.getdimensions())[1]
   k_cols = (k_gpu.getdimensions())[0]
   k_rows = (k_gpu.getdimensions())[1]
   gpuCopy, a_gpu, res_gpu
; padded arrays   
   gpuCopy, k_gpu, pad_gpu
   gpuAdd, 0.0, pad_gpu, 0.0, pad_gpu, 0.0, pad_gpu
   gpuSubArr, a_gpu, -1, -1, pad_gpu, [0,a_cols-1], [0,a_rows-1]  
   aHat_gpu = gpufft(pad_gpu,-1) 
   gpuCopy, k_gpu, pad_gpu
   fHat_gpu = gpufft(pad_gpu,-1)    
   gpuFree, pad_gpu   
   aHat_gpu = gpuMult(aHat_gpu,fHat_gpu,LHS=aHat_gpu)
   gpuFree, fHat_gpu
   aHat_gpu = gpufft(aHat_gpu,1,LHS=aHat_gpu) ; overwrites aHat_gpu
   aHatr_gpu = gpuReal(aHat_gpu)
   gpuFree, aHat_gpu
   gpuReform, aHatr_gpu, k_cols, k_rows
   gpuSubArr, aHatr_gpu, [0,a_cols-1], [0,a_rows-1], res_gpu, -1, -1
   gpuFree, aHatr_gpu   
   return, res_gpu
end