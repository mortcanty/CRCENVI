; docformat = 'rst'
; gpuem.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details
;+
; :Description:
;       Clustering with Gaussian mixtures, EM algorithm::
;           Gath, I. and Geva, A. B. (1989). Unsupervised 
;           optimal fuzzy clustering. 
;           IEEE Transactions on Pattern Analysis and 
;           Machine Intelligence, 3(3), 773–781.
; :Params:
;        G: in, required
;           observations data array (NNxn )
;           (row vectors)
;        U: in, out, required
;           initial class probability membership array (nxK)
;           (column vectors)
;        Ms: out, required  (NNxK)
;           cluster means matrix 
;           (row vectors)
;        Ps: out, required (K)
;           cluster priors
;        Cs: out, required 
;           cluster covariance matrices (KxNNxNN)
; :Keywords::
;        unfrozen: in, optional
;           indices of the observations which
;           take part in the iteration (default all)
;        wnd: in, optional
;           window for displaying the log likelihood
;        maxinter: in, optional
;           maximum iterations
;        mininter: in, optional
;           minimum iterations
;        pdens: out, optional
;           class partition densities 
;        fhv: out, optional
;           class fuzzy hypervolumes 
;        T0: in, optional
;           initial annealing temperature (default 1.0)
;        verbose: in, optional
;           set to print output info to IDL log
;           (default 0)
;        beta: in, optional
;           spatial field parameter
;        num_cols, num_rows: in, required
;            if beta>0 image dimensions 
;        status: out, optional
;            0 if successful else 1           
; :Uses:
;      COYOTE::
;      gpuconvol2D::
;      GPULIB                     
; :Author:
;      Mort Canty (2013) 
;-
Pro GPUEM, G, U, Ms, Ps, Cs, unfrozen=unfrozen, wnd=wnd, maxiter=maxiter, miniter=miniter, $
        pdens=pdens, fhv=fhv, T0=T0, verbose=verbose,$
        num_cols=num_cols, num_rows=num_rows, beta=beta, status=status
        
COMPILE_OPT IDL2

;resolve_routine, 'gpupow', /compile_full_file ; workaround for bug in gpuinit

common sd, seed 
        
status=0
sz = size(G)
n = (size(U))[1]
K = (size(U))[2]
if sz[1] eq 1 then message, 'GPUEM: data must have >1 dimension'
NN = sz[1]
; check keywords
if n_elements(unfrozen) eq 0 then unfrozen = lindgen(n)
if n_elements(maxiter) eq 0 then maxiter = 500
if n_elements(miniter) eq 0 then miniter = 10
if n_elements(T0) eq 0 then T0 = 1.0
if n_elements(beta) eq 0 then beta = 0.0
if beta gt 0 then begin
   if (n_elements(num_cols) eq 0) or (n_elements(num_rows) eq 0) then $
      message, 'GPUEM: beta>0 but no image dimensions supplied'
   if num_cols*num_rows ne n then message, 'GPUEM: beta>0 but wrong image dimensions'
endif
if n_elements(wnd) eq 1 then begin
   wset,wnd
   wxsize=!D.x_size
   wysize=!D.y_size
   if wxsize gt 400 and wysize gt 400 then $
      position=[0.2,0.9-400.0/wysize,0.2+400.0/wxsize,0.9] $
   else position = [0.2,0.1,0.9,0.9]
endif
if n_elements(verbose) eq 0 then verbose = 0

if verbose then begin
   print,'EM clustering '+strtrim(n_elements(unfrozen),2)+' pixels'
   print, 'Initial annealing temperature '+strtrim(T0)
endif

; covariance matrix
Cs = fltarr(NN,NN,K)
; log likelihood array
LL = fltarr(maxiter)
; partition density and fuzzy hypervolume arrays
pdens = (fhv = fltarr(K))

; GPU arrays
;  work
onesn_gpu = gpuPutArr(fltarr(n)+1)
onesNN_gpu = gpuPutArr(fltArr(NN)+1)
;  data matrix
G_gpu = gpuPutArr(float(G))
;  membership matrix
U_gpu = gpuPutArr(float(U))
unfrozen_gpu = gpuputArr(float(unfrozen))
;  filter kernel (padded and centered)
Nb = fltarr(num_cols+2,num_rows+2)
Nb[0:2,0:2] = [[0.0,0.25,0.0],[0.25,0.0,0.25],[0.0,0.25,0.0]]
Nb = shift(Nb,-1,-1)
Nb_gpu = gpuPutArr(Nb)

; --------- iteration ------------
progressbar = Obj_New('cgprogressbar', /cancel,title= 'EM iterating')
progressbar->start
dU = 1.0
iter=0L
T = T0
start_time = systime(2)
print, 'running CUDA ...'
while ((dU gt 0.001) or (iter lt miniter)) and (iter lt maxiter) do begin
   if progressbar->CheckCancel() then begin
      print,'clustering aborted'
      goto, done
   endif
   progressbar->Update,(iter*100)/maxiter
   
   gpuCopy, U_gpu, Uold_gpu
      
; number of pixels in each cluster
   ns_gpu = gpuTotal(U_gpu,1) 
   
; prior probabilities
   Ps_gpu = gpuAdd(1.0/n,ns_gpu,0.0,ns_gpu,0.0)
   
; cluster means   
   Ms_gpu = U_gpu##G_gpu
   tmp_gpu = gpuMatrix_Multiply(onesNN_gpu,ns_gpu,/btranspose)
   Ms_gpu = Ms_gpu/tmp_gpu   
   ns = gpuGetArr(ns_gpu)
   Ps = gpuGetArr(Ps_gpu)
   gpuFree, [tmp_gpu,ns_gpu,Ps_gpu]
   
; loop over the cluster index
   for j=0,K-1 do begin         
      gpuView,Ms_gpu,j*NN,NN,Msj_gpu
      gpuView,U_gpu,j*n,n,Uj_gpu
      
;    distances to jth mean      
      W_gpu = gpuMatrix_Multiply(Msj_gpu,onesn_gpu,/btranspose)
      Ds_gpu = G_gpu - W_gpu 
      
;    jth covariance matrix      
      Wt_gpu  = gpuTranspose(W_gpu)
      gpuFree,W_gpu 
      Dst_gpu = gpuTranspose(Ds_gpu)
      gpuFree,Ds_gpu   
      tmp_gpu = gpuSqrt(Uj_gpu)   
      for i=0,NN-1 do begin
         gpuView,Wt_gpu,i*n,n,Wti_gpu
         gpuView,Dst_gpu,i*n,n,Dsti_gpu    
         Wti_gpu = tmp_gpu*Dsti_gpu       
      endfor   
      gpuFree,tmp_gpu 
      C_gpu = gpuMatrix_Multiply(Wt_gpu,Wt_gpu,/atranspose) 
      gpuFree,Wt_gpu
      C_gpu = gpuAdd(1.0/ns[j],C_gpu,0.0,C_gpu,0.0)
      
 ;   invert it on the CPU in double precision   
      C = gpuGetArr(C_gpu)     
      gpuFree, C_gpu
      Cs[*,*,j]=C     
      sqrtdetC = sqrt(determ(C,/double))
      Cinv = invert(C,/double)
      gpuPutArr,float(Cinv),Cinv_gpu   
         
;    new memberships     
      gpuMatrix_Multiply,Dst_gpu,Cinv_gpu,tmp_gpu
      gpuMult,Dst_gpu,tmp_gpu,tmp_gpu
      qf_gpu = gpuTotal(tmp_gpu,2) 
      gpuFree,Cinv_gpu
      gpuFree,Dst_gpu
      gpuFree,tmp_gpu
      qfu_gpu = gpuSubscript(unfrozen_gpu,qf_gpu)
      gpuExp,Ps[j]/sqrtdetC,-0.5,qfu_gpu,0.0,0.0,tmp_gpu
      gpuFree,qfu_gpu
      gpuSubscript,unfrozen_gpu,tmp_gpu,Uj_gpu,/LEFTSUBSCRIPT
      gpuFree,tmp_gpu
         
;    class hypervolume and partition density
      fhv[j] = sqrtdetC     
      gpuLt,qf_gpu,onesn_gpu,qflt_gpu    
      gpuFree,qf_gpu 
      indices_gpu = gpuWhere(qflt_gpu,count)
      gpuFree,qflt_gpu
      if count gt 0 then begin
         Ujind_gpu = gpuSubscript(indices_gpu,Uj_gpu)
         pdens[j] = gpuTotal(Ujind_gpu)/fhv[j]
         gpuFree,Ujind_gpu
      endif
      gpuFree,indices_gpu
      
;    random membership for annealing
      if T gt 0.0 then begin
         gpuRandomu,seed,n,Ur_gpu
         pow_gpu = gpuMake_Array(n,/float,value=1.0/T)
         Ur_gpu = gpuPow(Ur_gpu,pow_gpu,LHS=Ur_gpu)
         gpuAdd, -1.0, Ur_gpu, 0.0, Ur_gpu,1.0, Ur_gpu
         gpuMult, Uj_gpu, Ur_gpu, Ur_gpu
         Ursub_gpu = gpuSubscript(unfrozen_gpu, Ur_gpu)
         gpuSubscript, unfrozen_gpu, Ursub_gpu, Uj_gpu, /LEFTSUBSCRIPT
         gpuFree, Ur_gpu
         gpuFree, pow_gpu
         gpuFree, Ursub_gpu 
      endif
   endfor 
      
; determine spatial membership matrix
   if beta gt 0 then begin
      for j=0,K-1 do begin
;    spatial membership        
;         U_N=1.0-convol(reform(U[*,j],num_cols,num_rows),Nb)
;         V[*,j] = exp(-beta*U_N)  
         gpuView, U_gpu, j*n, n, Uj_gpu
         gpuReform, Uj_gpu, num_cols, num_rows
         Vj_gpu = gpuConvol2D(Uj_gpu, Nb_gpu)
         gpuAdd, beta, Vj_gpu, 0.0, Vj_gpu, -beta, Vj_gpu
         gpuExp, Vj_gpu, Vj_gpu 
;    combine spectral/spatial for unfrozen
;         U[unfrozen,j]=U[unfrozen,j]*V[unfrozen,j]            
         gpuMult, Vj_gpu, Uj_gpu, Vj_gpu
         Vjunf_gpu = gpuSubscript(unfrozen_gpu, Vj_gpu)
         gpuSubscript, unfrozen_gpu, Vjunf_gpu,Uj_gpu, /LEFTSUBSCRIPT   
         gpuFree, Vj_gpu     
      endfor      
      gpuFree, Vjunf_gpu
   endif     
; normalize all   
   gpuAdd, 1.0, U_gpu, 0.0, U_gpu, 1e-20, U_gpu
   a_gpu = gpuTotal(U_gpu,2)
   for j=0,K-1 do begin
      gpuView, U_gpu, j*n, n, Uj_gpu
      gpuDiv, Uj_gpu, a_gpu, Uj_gpu
   endfor
   gpuFree,a_gpu   
                                    
; log likelihood
   gpuLog, U_gpu, logU_gpu
   gpuMult, Uold_gpu, logU_gpu, logU_gpu
   LL[iter] = gpuTotal(logU_gpu)  
   gpuFree, logU_gpu
   
; stopping criterion  
   gpuSub, U_gpu, Uold_gpu, Uold_gpu
   gpuAbs, Uold_gpu, Uold_gpu
   du = gpuMax(Uold_gpu)  
   gpuFree, Uold_gpu
      
   iter++
   T = 0.8*T
   if (n_elements(wnd) eq 1) and (iter mod 5 eq 0) then $
     plot, LL[0:iter-1], color=0, background= 'FFFFFF'XL, xtitle='Iteration',position=position
          
endwhile
; ---------- end iteration -----------
done:
if verbose then begin
   if iter lt maxiter then print, 'converged or interrupted after '+strtrim(iter,2)+' iterations' $
      else print, 'no convergence after '+strtrim(maxiter,2)+' iterations'
   print,'partition density for '+strtrim(K,2)+' clusters: ', total(pdens*fhv)/total(fhv)
   print, 'elapsed time ',systime(2)-start_time
   print, 'time per iteration ',(systime(2)-start_time)/iter
endif
gpuGetArr, U_gpu, U
gpuGetArr, Ms_gpu, Ms
; cleanup
gpuFree,[U_gpu,unfrozen_gpu,Nb_gpu,G_gpu,Ms_gpu,onesn_gpu,onesNN_gpu]
progressbar->destroy

End
