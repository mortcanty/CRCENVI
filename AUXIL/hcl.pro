; docformat = 'rst'
; hcl.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details

function best_merge, D
; return the best merge, larger index first
   n = (size(D))[1]
   _ = min(D,label)
   return, [label[0] mod n,label[0]/n]
end

;+
; :Description:
;       Agglomerative hierarchic clustering with
;       sum of squares cost function::
;          Fraley, C. (1996). Algorithms for model-based 
;          Gaussian hierarchical clustering.
;          Technical report 311, Department of Statistics, 
;          University of Washington, Seattle.
; :Params:
;        G: in, required  
;           data matrix
;        KK: in, required  
;            number of clusters
;        Ls: out, required 
;           cluster labels of observations  
; :Uses:
;      COYOTE                     
; :Author:
;      Mort Canty (2009) 
;-
pro hcl, G, KK, Ls

m = (size(G))[2]        ; number of sample vectors
mults = lonarr(m)+1     ; initial multiplicities
Ls = findgen(m)         ; initial cluster labels

; initialize merge-cost array in upper triangle
Delta = transpose(total(G^2,1))##(intarr(m)+1)
Delta = Delta + transpose(intarr(m)+1)##total(G^2,1)
; overwrite diagonal and lower triangle with 10e10 
Delta = Delta - 2*G##transpose(G)
idx = lindgen(m^2)
idx = where(idx mod m le idx/m)
Delta[idx] = 10e10

; begin iteration
cost = 0.0
c = m
progressbar = Obj_New('progressbar', Color='blue', Text='0',$
              title='HCL: classes remaining...',xsize=250,ysize=20)
progressbar->start
while c gt KK do begin
   if progressbar->CheckCancel() then begin
      print,'clustering aborted'
      Ls=-1
      progressbar->Destroy
      return
   endif
   progressbar->Update,(m-c)*100/m,text=strtrim(c,2)
   bm = best_merge(Delta)
   j = bm[0]
   i = bm[1]
; j > i   
   cost = cost + Delta[j,i]
; re-label   
   idx = where(Ls eq j,count)
   if count gt 0 then Ls[idx] = i   
   idx = where(Ls gt j,count)
   if count gt 0 then Ls[idx] = Ls[idx]-1
; pre-merge multiplicities
   ni = mults[i]
   nj = mults[j] 
; update merge-cost array, k = i+1 ... c-1 
   if c-i-1 eq 0 then k = [i+1] else k = i+1+lindgen(c-i-1)  
   nk = mults[k]
   Djk = Delta[j,k]<Delta[k,j]  
   idx = where(k eq j,count) 
   if count gt 0 then Djk[idx]=0 
   Delta[k,i] = ( (ni+nk)*Delta[k,i]+(nj+nk)*Djk-nk*Delta[j,i] )/(ni+nj+nk)  
; update merge-cost array, k = 0 ... i-1   
   if i eq 0 then k = [0] else k = lindgen(i)
   nk = mults[k]
   Djk = Delta[j,k]<Delta[k,j]
   idx = where(k eq j,count) 
   if count gt 0 then Djk[idx]=0  
   Delta[i,k] = ( (ni+nk)*Delta[i,k]+(nj+nk)*Djk-nk*Delta[j,i] )/(ni+nj+nk)  
; update multiplicities
   mults[i] = mults[i]+mults[j]
; delete the upper cluster   
   idx = replicate(1L,c)
   idx[j]=0L
   idx = where(idx eq 1L)
   mults = mults[idx]
   Delta = Delta[idx,*]
   Delta = Delta[*,idx]
   c = c-1
endwhile
progressbar->destroy
end
