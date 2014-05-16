; docformat = 'rst'
; gpukernelproject.pro
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
;      Centers and projects test data for
;      kernel PCA or kernel MAF (linear or gaussian)::     
;      returns 1 on success else 0        
; :Params:
;      alpha: in, required 
;         normalized eigenvectors of kernel matrix
;      gma: in, required 
;         radial basis kernel parameter
;      G: in, required  
;         training data matrix
;      GG: in, required  
;         test data matrix to be projected
;      num_cols, num_rows: in, required 
;         image dimensions
;      image: out
;         projected data as BSQ image  
; :Keywords:
;      center_train: in, optional
;         if set, center on training means (default 0) 
;      nscale: in, optional,type=float
;            multiple of scale for GMA when KERNEL=1 (default 1.0)    
;      kernel: in, optional, type=integer
;           the kernel used
;            0: linear
;            1: Gaussian (default)         
; :Uses:
;      ENVI::
;      CENTER::
;      GPUKERNEL_MATRIX::    
;      COYOTE::
;      GPULIB
; :Author:
;      Mort Canty (2009)       
;-
function gpukernelProject, alpha, gma, G, GG, num_cols, num_rows, image, $
   kernel=kernel, center_train=center_train, nscale=nscale

COMPILE_OPT STRICTARR

; kernel matrix blocks for row-wise processing 
;   <--------------- num_pixels ------------------->   -
;  (            |            |  ...  |              )  |
;  (            |            |  ...  |              )  |
;  (    block   |   block    |  ...  |    block     )  m
;  (      0     |     1      |  ...  |  num_rows-1  )  |
;  (            |            |  ...  |              )  |
;  (            |            |  ...  |              )  -
;                <-num_cols->

if n_elements(center_train) eq 0 then center_train=0 else center_train=center_train
if n_elements(kernel) eq 0 then kernel = 1 else kernel=kernel
if n_elements(nscale) eq 0 then nscale = 1.0 else nscale=nscale

m = n_elements(G[0,*])
num_bands = n_elements(G[*,0])
num_comps = n_elements(alpha[0,*])
num_pixels = n_elements(GG[0,*])
progressbar = Obj_New('cgprogressbar', /cancel, $
             title='Kernel transform ',xsize=300,ysize=20)
progressbar->start
if gpu_detect() then begin
; using CUDA: put dual variables (normalized eigenvectors) onto GPU 
   alpha_gpu = gpuPutArr(alpha)
; put other arrays onto GPU 
   image_gpu = gpuFltArr(num_cols*num_comps,num_rows)  ;first num_comps projections, row-by-row  
   G_gpu = gpuPutArr(G)                                ;training data matrix
   GG_gpu = gpuPutArr(GG)                              ;data matrix
   result=1
   if not center_train then begin
; centering on test means   
      rowmeans_gpu = gpuFltArr(m)                         ;work arrays  
      onesnc_gpu = gpuputarr(fltarr(num_cols)+1.0/num_pixels)  
;    centering
      meanKK = 0.0       
      i = 0L 
      while i lt num_rows do begin
         pct=i*100/num_rows
         progressbar->Update,fix(pct)
;       ggi = GG[*,i*num_cols:(i+1)*num_cols-1]
         gpuView,GG_gpu,i*num_cols*num_bands,num_bands*num_cols,GGi_gpu       
         gpuReform,GGi_gpu,num_bands,num_cols 
         KK_gpu = gpukernel_matrix(G_gpu,GGi_gpu,kernel=kernel,gma=gma,nscale=nscale)
         meanKK = meanKK + gputotal(KK_gpu)/(m*num_pixels)
;       row means
         tmp_gpu = gpumatrix_multiply(onesnc_gpu,KK_gpu,/atranspose)
         rowmeans_gpu = gpuadd(rowmeans_gpu,tmp_gpu,lhs=rowmeans_gpu)
         gpufree,KK_gpu
         gpuFree,tmp_gpu 
         if progressbar->CheckCancel() then begin
             print,'aborted...'
             i=num_rows-1
             result=0
         endif
         i++ ; next row
      endwhile
      gpuFree, onesnc_gpu
      txt = 'Pass 2 of 2: projecting... '
   end else begin
   ; centering on training means  
      K_gpu = gpukernel_matrix(G_gpu,kernel=kernel,gma=gma,nscale=nscale)
      meanKK = gpuTotal(K_gpu)/(m*m)
;    row means
      rowmeans_gpu = gpuTotal(K_gpu,1)  
      gpuFree,K_gpu
      rowmeans_gpu = gpuAdd(1./m, rowmeans_gpu, 0.0, rowmeans_gpu, 0.0, LHS=rowmeans_gpu)   
      txt = 'Projecting... '
   endelse   
   onesm_gpu = gpuputarr(fltarr(m)+1.0)
   onesnc_gpu = gpuputarr(fltarr(num_cols)+1.0)     
   progressbar->destroy
   progressbar = Obj_New('cgprogressbar', /cancel,$
              title='Kernel transform ',xsize=300,ysize=20) 
   progressbar->start           
; projecting
   i = 0L
   while i lt num_rows do begin
      pct=i*100/num_rows
      progressbar->Update,fix(pct)
;    ggi = GG[*,i*num_cols:(i+1)*num_cols-1]
      gpuView,GG_gpu,i*num_cols*num_bands,num_bands*num_cols,GGi_gpu       
      gpuReform,GGi_gpu,num_bands,num_cols 
      KK_gpu = gpukernel_matrix(G_gpu,GGi_gpu,kernel=kernel,gma=gma,nscale=nscale)
;    subtract row means      
      tmp1_gpu = gpumatrix_multiply(onesnc_gpu,rowmeans_gpu,/btranspose)
      KK_gpu = gpusub(KK_gpu,tmp1_gpu,lhs=KK_gpu)
;    subtract column means     
      colmeans_gpu = gpuTotal(KK_gpu,2)
      colmeans_gpu = gpuAdd(1./m, colmeans_gpu, 0.0, colmeans_gpu, 0.0, LHS=colmeans_gpu) 
      tmp2_gpu = gpumatrix_multiply(colmeans_gpu,onesm_gpu,/btranspose)    
      gpuFree,colmeans_gpu  
      KK_gpu = gpusub(KK_gpu,tmp2_gpu,lhs=KK_gpu)
      gpufree,tmp1_gpu
      gpuFree,tmp2_gpu
;    add overall mean
      KK_gpu = gpuadd(1.0,KK_gpu,0.0,KK_gpu,meanKK,lhs=KK_gpu)        
;    image[*,i,*] = alpha##KK
      gpuView,image_gpu,i*num_cols*num_comps,num_cols*num_comps,rowi_gpu
      gpuReform,rowi_gpu,num_cols,num_comps
      rowi_gpu = gpuMatrix_multiply(KK_gpu,alpha_gpu,lhs=rowi_gpu) 
      gpufree,KK_gpu   
      if progressbar->CheckCancel() then begin
          print,'interrupted...'
          i=num_rows-1
          result=0
      endif
      i++ ; next row
   endwhile
; get result back to CPU
   image = transpose(reform(gpugetarr(image_gpu),num_cols,num_comps,num_rows), [0, 2, 1])
; tidy up   
   gpuFree,[image_gpu,G_gpu,GG_gpu,alpha_gpu,rowmeans_gpu,onesm_gpu,onesnc_gpu]  
   progressbar->destroy 
   return, result
end else begin
; CUDA not available
   image = fltarr(num_cols,num_rows,num_comps)  
   result=1 
   if not center_train then begin
; centering  
      rowmeans = fltarr(m)
      meanKK = 0.0 
      i = 0L  
      while i lt num_rows do begin
         pct=i*100/num_rows
         progressbar->Update,fix(pct)
;       get the ith row of the image
         GGi = GG[*,i*num_cols:(i+1)*num_cols-1]
;       m x num_cols matrix of kernels
         KK = kernel_matrix(G,GGi,kernel=kernel,gma=gma,nscale=nscale) 
         rowmeans = rowmeans + total(KK,1)/num_pixels
         meanKK = meanKK+total(KK)/(m*num_pixels)
         if progressbar->CheckCancel() then begin
            print,'aborted...'
            i=num_rows-1
            result=0
         endif
         i++ ; next row
      endwhile
      txt = 'Pass 2 of 2: projecting... '
   end else begin  
;    centering on training means   
      K = kernel_matrix(G,kernel=kernel,gma=gma,nscale=nscale)
      rowmeans = total(K,1)/m
      meanKK = total(K)/(m*m)
      txt = 'Projecting... '
   endelse 
; projecting   
   i = 0L
   progressbar->destroy
   progressbar = Obj_New('cgprogressbar', /cancel,$
              title='Kernel transform ',xsize=300,ysize=20) 
   progressbar->start      
   while i lt num_rows do begin
      pct=i*100/num_rows
      progressbar->Update,fix(pct)
;    get the ith row of the image
      GGi = GG[*,i*num_cols:(i+1)*num_cols-1]
;    m x num_cols matrix of kernels
      KK = kernel_matrix(G,GGi,kernel=kernel,gma=gma,nscale=nscale)
;    subtract column means   
      colmeans = total(KK,2)/m   
      KK = KK - (fltarr(1,m)+1)##colmeans
;    subtract row means
      KK = KK - transpose(rowmeans)##(fltarr(num_cols)+1)
;    add overall mean
      KK = KK + meanKK          
;    project
      image[*,i,*] = alpha##KK
      if progressbar->CheckCancel() then begin
         print,'interrupted...'
         i=num_rows-1
         result=0
      endif
      i++ ; next row
   endwhile   
   progressbar->destroy
   return, result
endelse   

end