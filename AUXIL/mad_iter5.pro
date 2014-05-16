; docformat = 'rst'
; mad_iter5.pro
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
;       Function for Iteratively Re-weighted Multivariate
;       Alteration Detection (IR-MAD)::
;          Nielsen, A. A. (2007). The regularized 
;          iteratively reweighted MAD method for change 
;          detection in multi- and hyperspectral data.
;          IEEE Transactions on Image Processing, 
;          16(2), 463–478.
; :Params:
;       fid1,fid2: in,required
;          image file IDs
;       dims1,dims2: in,required
;          image spatial subsets
;       pos1,pos2: in,required
;          image spectral subsets
;       m_fid: in,required (-1 for no masking)   
;          mask band ID      
; :KEYWORDS:
;       A, B: out, optional 
;          arrays of transformation eigenvectors in columns,
;          smallest correlation first
;       means1,means2: out,optional   
;          weighted mean values, row-replicated
;       rho: out,optional
;         canonical correlations in decreasing order
;       sigma: out,optional
;         standard deviations of MAD variates
;         in increasing order
;       lam: in,optional
;            length penalization parameter
;            (default 0.0)
;       niter: in,optional
;            max number of iterations,
;            (default is 50)
;       verbose: in,optional
;           print intermediate results to log
;           (default 0);       
; :Uses:
;       ENVI:: 
;       CPM_DEFINE::      
;       GEN_EIGENPROBLEM::       
;       COYOTE
; :Author:
;       Mort Canty (2012)        
;-
Function mad_iter5, fid1, fid2, dims1, dims2, pos1, pos2, m_fid, niter=niter, rho=rho, $
               sigma=sigma, A=A, B=B, means1=means1, means2=means2, $
               verbose=verbose, lam=lam

   COMPILE_OPT idl2                

   if n_elements(niter) eq 0 then niter = 50
   if n_elements(verbose) eq 0 then verbose = 0
   if n_elements(lam) eq 0 then lam = 0.0

   envi_file_query, fid1, interleave=interleave1
   envi_file_query, fid2, interleave=interleave2

   num_cols = dims1[2]-dims1[1]+1
   num_rows = dims1[4]-dims1[3]+1
   num_pixels = (num_cols*num_rows)
   num_bands = n_elements(pos1)
   
   if num_bands eq 1 then interleave1 = (interleave2 = 2)
   
; column vector of 1s   
   ones = fltarr(num_cols,1)+1.0

   if m_fid ne -1 then mask = envi_get_data(fid=m_fid,dims=dims1,pos=0) $
      else mask = bytarr(num_cols,num_rows)+1B
   if niter gt 0 then rhos=fltarr(num_bands,niter)-1

   if niter gt 1 then begin
      window, 12, xsize=600, ysize=400, title='IR-MAD'
      wset,12
      plot, indgen(niter), rhos[0,*], psym=4, xrange = [0,niter+1], yrange = [0,1.1], $
         xtitle='iteration', ytitle='Correlations',color=0,background='FFFFFF'XL
   endif

; get tile ids
   tile_id1 = envi_init_tile(fid1,pos1,num_tiles=num_tiles, $
       interleave=interleave1,xs=dims1[1],xe=dims1[2],ys=dims1[3],ye=dims1[4])
   tile_id2 = envi_init_tile(fid2,pos2, $
       interleave=interleave2,xs=dims2[1],xe=dims2[2],ys=dims2[3],ye=dims2[4])

   sigMADs = fltarr(num_bands,num_cols)
   
; length penalization
   Omega_L = Identity(num_bands)

;; slope penalization
;   L = fltarr(num_bands,num_bands-1)
;   row = fltarr(num_bands)
;   row[0:1] = [1,-1]
;   for i=0,num_bands-2  do L[*,i] = shift(row,i)
;   Omega_S = transpose(L)##L
;
;; curvature penalization
;   L = fltarr(num_bands,num_bands-2)
;   row = fltarr(num_bands)
;   row[0:2] = [1,-2,1]
;   for i=0,num_bands-3  do L[*,i] = shift(row,i)
;   Omega_C = transpose(L)##L

; begin iteration
   iter=0L
   interrupting = 0
   delta = 1.0
   old_rho = fltarr(num_bands)
   repeat begin
      progressbar = Obj_New('cgprogressbar', Text=' ',$
                    title='IR-MAD',xsize=300,ysize=20,/cancel)
      progressbar->start
      cpm = obj_new('CPM',2*num_bands)
      txt = 'Iter '+strtrim(iter,2)
; spectral tiling
      for tile_index=0L,num_tiles-1 do begin
         if progressbar->CheckCancel() then begin
            if (iter gt 0) and not interrupting then begin
               iter = niter
               interrupting = 1
               txt = 'Interrupting...'
            end else begin
               print,'Calculation aborted'
               obj_destroy, cpm
               envi_tile_done, tile_id1
               envi_tile_done, tile_id2
               wdelete,12
               message, 'Calculation aborted'
            endelse
         endif
         if tile_index mod 10 eq 0 then begin
            pct = (tile_index)*100/(num_tiles)
            progressbar->Update,fix(pct)
         endif
         tile1 = envi_get_tile(tile_id1,tile_index)
         tile2 = envi_get_tile(tile_id2,tile_index)
         if interleave1 eq 1 then tile1 = transpose(tile1)
         if interleave2 eq 1 then tile2 = transpose(tile2)
         if iter gt 0 then begin
            mads = ((tile1-means1)##A-(tile2-means2)##B)
            chi_sqr = total((mads/sigMADs)^2,1)
            weights=1-chisqr_pdf(chi_sqr,num_bands)
         end else weights = ones   
;       only sample pixels under the mask
         indices = where(mask[*,tile_index],count)
         if count gt 0 then $
            cpm->update,[tile1[*,indices],tile2[*,indices]],weights=weights[indices]
      endfor
      progressbar->destroy

; canonical correlation
      SS = cpm->Covariance()
      means = cpm->Means()
      Obj_Destroy, cpm
     
      S11 = SS[0:num_bands-1,0:num_bands-1]
      S11 = (1-lam)*S11 + lam*omega_L
      S22 = SS[num_bands:*,num_bands:*] 
      S22 = (1-lam)*S22 + lam*omega_L
      S12 = SS[num_bands:*,0:num_bands-1]
      S21 = SS[0:num_bands-1,num_bands:*]
 
      C1 = S12 ## invert(S22,/double) ## S21
      B1 = S11
      C2 = S21 ## invert(S11,/double) ## S12
      B2 = S22
      
      if num_bands gt 1 then begin
         gen_eigenproblem, C1, B1, A, mu2 
         gen_eigenproblem, C2, B2, B, mu2
      end else begin
         mu2 = [C1/B1]
         A = [1/sqrt(B1)]
         B = [1/sqrt(B2)]
      endelse         
      
      mu = sqrt(mu2)  
      a2=diag_matrix(transpose(A)##A)
      b2=diag_matrix(transpose(B)##B) 
      sigma = sqrt( (2-lam*(a2+b2))/(1-lam)-2*mu )
      rho=mu*(1-lam)/sqrt( (1-lam*a2)*(1-lam*b2) )     
           
      sigMads = ones##sigma 
      means1  = ones##means[0:num_bands-1]
      means2  = ones##means[num_bands:*]
      
; stopping criterion
      delta = max(abs(rho-old_rho))
      old_rho = rho

; ensure sum of positive correlations between X and U is positive
; their covariance matrix is S11##A
      invsig_x = diag_matrix(1/sqrt(diag_matrix(S11)))
      sum = total(invsig_x##S11##A,2)
      A = A##diag_matrix(sum/abs(sum))   

; ensure positive correlation between each pair of canonical variates
      cov = diag_matrix(transpose(A)##S12##B)
      B = B##diag_matrix(cov/abs(cov))

      if iter gt 0 and iter eq niter then goto, done

;    print sigmas and plot canonical correlations
      if verbose then begin
         print, 'delta = '+strtrim(delta,2)
         print, 'iteration '+strtrim(iter,2)
         print, reverse(sigMADs[*,0]), format='("Sigma MADs",'+strtrim(num_bands,2)+'f10.5)'
      endif
      if (niter gt 1) and (iter lt niter) then begin
         rhos[*,iter] = rho
         wset,12
         plot, indgen(iter+1)+1, rhos[0,0:iter], psym=-4, xrange=[0,niter+1],  yrange = [0,1.1], $
           xtitle='iteration', ytitle='Correlations',color=0,background='FFFFFF'XL
         for i=1,num_bands-1 do $
           oplot, indgen(iter+1)+1, rhos[i,0:iter], psym=-4, color=0
      endif
   done:
      iter=iter+1
   endrep until (iter gt niter) or (delta lt 0.001)
   envi_tile_done, tile_id1
   envi_tile_done, tile_id2

; successfully finished, so
   return, 0

End