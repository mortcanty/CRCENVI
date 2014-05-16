; docformat = 'rst'
; gpuCov__define.pro
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
;       Object class for variance/covariance matrix calculation on GPU      
; :Params:
;       NN, n:  in, required, type=longint 
;          dimensions of the data matrices                              
; :Uses:
;       GPULib
; :Author:
;       Mort Canty (2009) 
;-
function GPUCOV::Init, NN,n
   self.n = n
   self.NN = NN
   self.one_n_gpu = ptr_new(gpuPutArr(fltArr(n)+1)) 
   self.one_NN_gpu = ptr_new(gpuPutArr(fltArr(NN)+1))
   return,1
end   

pro GPUCOV::cleanup
   gpuFree, [*self.one_n_gpu,*self.one_NN_gpu]
   ptr_free, self.one_n_gpu
   ptr_free, self.one_NN_gpu
end

;+ 
; :Params:
;       G_gpu:  in, required, type={GPUHANDLE} 
;          handle to a NNxn-dimensional data matrix    
; :Keywords:
;       weights_gpu: in, optional, type ={GPUHANDLE}
;          handle to n-dimensional array of weights  
;-
function GPUCOV::Covariance, G_gpu, weights_gpu=ws_gpu 
   if n_elements(ws_gpu) eq 0 then begin
      Ms_gpu = gpuTotal(G_gpu,2)
      gpuAdd, 1.0/self.n, Ms_gpu, 0.0, Ms_gpu, 0.0, Ms_gpu
      gpuMatrix_Multiply, Ms_gpu, *self.one_n_gpu, ds_gpu, /btranspose
      gpuSub,G_gpu,ds_gpu,ds_gpu
      gpuMatrix_Multiply,ds_gpu,ds_gpu,result_gpu,/btranspose 
      gpuAdd, 1.0/(self.n-1),result_gpu,0.0,result_gpu,0.0,result_gpu
      gpuFree,[Ms_gpu,ds_gpu]
   end else begin
 ; ensure that weights are normalized
      sum = gpuTotal(ws_gpu)
      gpuAdd, 1.0/sum, ws_gpu, 0.0, ws_gpu, 0.0, ws_gpu   
 ; weighted means   
      gpuMatrix_Multiply,*self.one_NN_gpu,ws_gpu,tmp0_gpu,/bTranspose 
      gpuMult,G_gpu,tmp0_gpu,tmp1_gpu
      Ms_gpu = gpuTotal(tmp1_gpu,2) 
; weighted covariance mtrix      
      gpuMatrix_Multiply,Ms_gpu,*self.one_n_gpu,ds_gpu,/btranspose
      gpuSub,G_gpu,ds_gpu,ds_gpu 
      gpuSqrt,tmp0_gpu,tmp0_gpu
      gpuMult,tmp0_gpu,ds_gpu,ds_gpu
      gpuMatrix_Multiply,ds_gpu,ds_gpu,result_gpu,/btranspose
      gpuAdd,sum/(sum-1),result_gpu,0.0,result_gpu,0.0,result_gpu
      gpuFree, [Ms_gpu,tmp0_gpu,tmp1_gpu,ds_gpu]
   endelse
   return,result_gpu
end   

Pro gpuCov__Define
class = {GPUCOV, $
         n:  0L, $
         NN: 0L, $
         one_n_gpu:  ptr_new(), $
         one_NN_gpu: ptr_new()}
End         