; docformat = 'rst'
; fkm.pro
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
;       Fuzzy Kmeans clustering algorithm::
;          Dunn, J. C. (1973). A fuzzy relative of the
;          isodata process and its use in detecting 
;          compact well-separated clusters. 
;          Journal of Cybernetics, PAM1-1, 32–57.
; :Params:
;        G: in, required
;           input data matrix
;        U: in, out, required
;           class probability membership matrix
;        Ms: out, required
;           cluster means (output) 
; :Uses:
;      COYOTE                     
; :Author:
;      Mort Canty (2009) 
;-
pro FKM, G, U, Ms, niter=niter, unfrozen=unfrozen

n = (size(G))[2]
NN = (size(G))[1]
K = (size(U))[2]

if n_elements(niter) eq 0 then niter = 500
if n_elements(unfrozen) eq 0 then unfrozen = lindgen(n)

; vector distances to cluster centers
Ds = fltarr(n,NN)
; work array
W = fltarr(n,NN)

; iteration
progressbar = Obj_New('progressbar', Color='blue', Text='0',$
              title='FKM clustering: delta_U...',xsize=250,ysize=20)
progressbar->start
dU = 1.0 & iter=0L
while ((dU gt 0.001) or (iter lt 20)) and (iter lt niter) do begin
   if progressbar->CheckCancel() then begin
      print,'clustering aborted'
      progressbar->Destroy
      return
   endif
   progressbar->Update,(iter*100)/niter,text=strtrim(dU,2)
   Uold = U
   UU = U*U
   ; update means and distances
   Ms = transpose(UU##G)
   for j=0,K-1 do begin
      Ms[j,*]=Ms[j,*]/total(UU[*,j])
      for i=0,NN-1 do W[*,i]=replicate(Ms[j,i],n)
      Ds = transpose(G)-W
      dd = 1/total(Ds*Ds,2)
      U[unfrozen,j] = dd[unfrozen]
   endfor
; normalize
   sums = total(U,2)
   for j=0,K-1 do U[*,j]=U[*,j]/sums
   dU = max(abs(U-Uold))
   iter=iter+1
endwhile

Ms = transpose(Ms)
progressbar->destroy

end
