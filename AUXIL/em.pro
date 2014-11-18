; docformat = 'rst'
; em.pro
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
;           observations data matrix
;        U: in, out, required
;           initial class probability membership matrix
;           (column vectors)
;        Ms: out, required
;           cluster means 
;        Ps: out, required
;           cluster priors
;        Cs: out, required
;           cluster covariance matrices 
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
;           partition density 
;        fhv: out, optional
;           fuzzy hypervolume 
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
;      INTERSECTION::
;      COYOTE                     
; :Author:
;      Mort Canty (2013) 
;-
Pro EM, G, U, Ms, Ps, Cs, unfrozen=unfrozen, wnd=wnd, maxiter=maxiter, miniter=miniter, $
        pdens=pdens, fhv=fhv, T0=T0, verbose=verbose, $
        num_cols=num_cols, num_rows=num_rows, beta=beta, status=status
common sd, seed  

COMPILE_OPT IDL2

catch, theError 
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
   progressbar->destroy
   status=1
   return
endif             
        
status=0
V = U*0.0
sz = size(G)
n = (size(U))[1]
qf = fltarr(n)
K = (size(U))[2]
if sz[0] eq 1 then NN = 1 else NN = sz[1]
Nb = [[0.0,0.25,0.0],[0.25,0.0,0.25],[0.0,0.25,0.0]]
; check keywords
if n_elements(unfrozen) eq 0 then unfrozen = lindgen(n)
if n_elements(maxiter) eq 0 then maxiter = 500
if n_elements(miniter) eq 0 then miniter = 10
if n_elements(T0) eq 0 then T0 = 1.0
if n_elements(beta) eq 0 then beta = 0.0
if beta gt 0 then begin
   if (n_elements(num_cols) eq 0) or (n_elements(num_rows) eq 0) then $
      message, 'EM: beta>0 but no image dimensions supplied'
   if num_cols*num_rows ne n then message, 'EM: beta>0 but wrong image dimensions'
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
if NN eq 1 then Cs = dblarr(K) else Cs = dblarr(NN,NN,K)
; vector distances to cluster centers
if NN eq 1 then Ds = dblarr(n) else Ds = dblarr(NN,n)
; work array
if NN eq 1 then W = dblarr(n) else W = dblarr(NN,n)
; log likelihood array
LL = fltarr(maxiter)
; partition density and fuzzy hypervolume arrays
pdens = (fhv = fltarr(K))

; iteration
progressbar = Obj_New('cgprogressbar',/cancel, $
              title= 'EM iterating')
progressbar->start
dU = 1.0
iter=0L
T = T0
start_time = systime(2)
while ((dU gt 0.001) or (iter lt miniter)) and (iter lt maxiter) do begin
  if progressbar->CheckCancel() then begin
     print,'clustering aborted'
     goto, done
  endif
  progressbar->Update,(iter*100)/maxiter
  Uold = U
; number of pixels in each cluster
  ns = total(U,1)
; prior probabilities
  Ps = ns/n
; calculate cluster means
; note: A ## TRANSPOSE(B) is equivalent to MATRIX_MULTIPLY(B, A, /ATRANSPOSE)
  Ms = matrix_multiply(U, G, /atranspose,/btranspose)
; loop over the cluster index 
  if NN gt 1 then for j=0,K-1 do begin
    Ms[j,*] = temporary(Ms[j,*])/ns[j]
    W = (fltarr(n)+1)##transpose(Ms[j,*])
    Ds = G - W   
;  covariance matrix
    for i=0,NN-1 do W[i,*] = $
         sqrt(U[*,j])*Ds[i,*]
    C = matrix_multiply(W,W,/btranspose)/ns[j]
    Cs[*,*,j] = C
    sqrtdetC = sqrt(determ(C,/double))
    Cinv = invert(C,/double)
    qf = total(Ds*(Ds##Cinv),1)    
;  class hypervolume and partition density
    fhv[j] = sqrtdetC
    indices = where(qf lt 1.0, count)
    if (count gt 0) then begin
       pdens[j] = total(U[indices,j])/fhv[j]
    endif    
;  new memberships
    U[unfrozen,j] = $
          exp(-qf[unfrozen]/2.0)*(Ps[j]/sqrtdetC)         
;  random membership for annealing
    if T gt 0.0 then begin
      Ur=1.0-randomu(seed,n_elements(unfrozen))^(1.0/T)
      U[unfrozen,j] = temporary(U[unfrozen,j])*Ur
    endif
  endfor $
; loop over the cluster index for the case NN=1
  else for j=0,K-1 do begin
     Ms[j] = Ms[j]/ns[j]
     replicate_inplace, W, Ms[j]
     Ds = G-W
;   variance
     Cs[j] = total(U[*,j]*Ds^2)/ns[j]
     sqrtdetC = sqrt(Cs[j])
;   class hypervolume and partition density
     fhv[j] = sqrtdetC
     indices = where(Ds^2/Cs[j] lt 1.0, count)
     if count gt 0 then pdens[j] = total(U[indices,j])/fhv[j]
;   new memberships
     U[unfrozen,j] = Ps[j]*exp(-0.5*Ds[unfrozen]^2/Cs[j])/sqrtdetC
;   random membership for annealing
     if T gt 0.0 then begin
       Ur = 1.0 - randomu(seed,n_elements(unfrozen))^(1.0/T)
       U[unfrozen,j] = temporary(U[unfrozen,j])*Ur
     endif
  endfor
; determine spatial membership matrix
  if beta gt 0 then begin
    for j=0,K-1 do begin
      U_N=1.0-convol(reform(U[*,j],num_cols, $
                            num_rows),Nb,/center)
      V[*,j] = exp(-beta*U_N)
    endfor
;  combine spectral/spatial for unfrozen
    U[unfrozen,*]=temporary(U[unfrozen,*])* $
                            V[unfrozen,*]
  endif                                        
; normalize all
  a = 1/total(U,2)
  void=where(finite(a),complement=complement, $
    ncomplement=ncomplement)
  if ncomplement gt 0 then a[complement]=0.0
  for j=0,K-1 do U[*,j]=temporary(U[*,j])*a
  
; log likelihood
  indices=where(U,count)
  if count gt 0 then LL[iter] =  $
       total(Uold[indices]*alog(U[indices]))
  dU = max(abs(U-Uold))
  iter=iter+1
  T = 0.8*T
  if (n_elements(wnd) eq 1) and (iter mod 5 eq 0) then $
    plot, LL[0:iter-1], color=0, background= 'FFFFFF'XL, xtitle='Iteration',position=position
endwhile
done:
Ms = transpose(Ms)
if verbose then begin
   if iter lt maxiter then print, 'converged or interrupted after '+strtrim(iter,2)+' iterations' $
      else print, 'no convergence after '+strtrim(maxiter,2)+' iterations'
   print,'partition density for '+strtrim(K,2)+' clusters: ', total(pdens)/total(fhv)
   print, 'elapsed time ',systime(2)-start_time
   print, 'time per iteration ',(systime(2)-start_time)/iter
endif
progressbar->destroy

End
