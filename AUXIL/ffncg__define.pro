; docformat = 'rst'
; ffncg__define.pro
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
;     Object class for implementation of a
;     two-layer, feed-forward neural network classifier.
;       Implements scaled conjugate gradient training::
;           Bishop, C. M. (1995). Neural Networks for 
;           Pattern Recognition. Oxford ;University Press.    
; :Inherits:
;        FFN__DEFINE       
; :Params:
;        Gs: in, required 
;           array of observation column vectors     
;        Ls: in, required
;            array of class label column vectors
;            of form (0,0,1,0,0,...0)^T          
;        L: in, required
;           number of hidden neurons
; :Examples:   
;     ffn = Obj_New("FFNCG",Gs,Ls,L)
; :Uses:
;       COYOTE    
; :Author:      
;       Mort Canty (2009)           
;-
Function FFNCG::Init, Gs, Ls, L, epochs=epochs
COMPILE_OPT STRICTARR
   catch, theError
   if theError ne 0 then begin
      catch, /cancel
      ok = dialog_message(!Error_State.Msg + ' Returning...', /error)
      return, 0
   endif
; initialize the superclass
   if not self->FFN::Init(Gs, Ls, L) then return, 0
   if n_elements(epochs) eq 0 then self.epochs = 1000L $
      else self.epochs = epochs
   self.cost_array = ptr_new(fltarr(self.epochs))
   self.valid_cost_array = ptr_new(fltarr(self.epochs))
   return, 1
End

Pro FFNCG::Cleanup
   ptr_free, self.cost_array
   ptr_free, self.valid_cost_array
   self->FFN::Cleanup
End

Function FFNCG::Gradient
COMPILE_OPT STRICTARR
   M = self->vForwardPass(*self.GTs,N)
   D_o = *self.LTs - M
   D_h = (N*(1-N)*((*self.Wo)##D_o))[*,1:*]
   dEh = -(*self.GTs)##transpose(D_h)
   dEo = -N##transpose(D_o)
   return, [dEh[*],dEo[*]]
End

Function FFNCG::Rop, V
COMPILE_OPT STRICTARR
; reform V to dimensions of Wh and Wo and transpose
  VhT=transpose(reform(V[0:self.LL*(self.NN+1)-1], $
                  self.LL,self.NN+1))
  Vo=reform(V[self.LL*(self.NN+1):*],self.KK,self.LL+1)
  VoT = transpose(Vo)
; transpose the weights
  Wo  = *self.Wo
  WoT = transpose(Wo)
; vectorized forward pass
  M = self->vForwardPass(*self.GTs,N)
; evaluation of v^T.H
  Zeroes = fltarr(self.ntp)
  D_o=*self.LTs-M                          ;d^o
  RIh=VhT##(*self.GTs)                     ;Rv{I^h}
  RN=N*(1-N)*[[Zeroes],[RIh]]              ;Rv{n}
  RIo=WoT##RN + VoT##N                     ;Rv{I^o}
  Rd_o=-M*(1-M)*RIo                        ;Rv{d^o}
  Rd_h=N*(1-N)*((1-2*N)*[[Zeroes],[RIh]]*(Wo##D_o) $
         + Vo##D_o + Wo##Rd_o)
  Rd_h=Rd_h[*,1:*]                         ;Rv{d^h}
  REo=-N##transpose(Rd_o)-RN##transpose(D_o);Rv{dE/dWo}
  REh=-*self.GTs##transpose(Rd_h)          ;Rv{dE/dWh}
  return, [REh[*],REo[*]]                  ;v^T.H
End

Function FFNCG::Hessian
COMPILE_OPT STRICTARR
   nw = self.LL*(self.NN+1)+self.KK*(self.LL+1)
   v = diag_matrix(fltarr(nw)+1.0)
   H = fltarr(nw,nw)
   for i=0,nw-1 do H[*,i] = self -> Rop(v[*,i])
   return, H
End

Function FFNCG::Eigenvalues
COMPILE_OPT STRICTARR
   H = self->Hessian()
   H = (H+transpose(H))/2
   return, eigenql(H,/double)
End

Pro FFNCG::Train, key=key
COMPILE_OPT STRICTARR
   if n_elements(key) eq 0 then key=0
; if validation option then train only on odd examples
   if key then void = $
   where((lindgen(self.np) mod 2) eq 0, $
                      complement=indices) $
   else indices = lindgen(self.np)
   self.ntp = n_elements(indices)
   self.GTs = ptr_new((*self.Gs)[indices,*])
   self.LTs = ptr_new((*self.Ls)[indices,*])
   w = [(*self.Wh)[*],(*self.Wo)[*]]
   nw = n_elements(w)
   g = self->gradient()
   d = -g    ; search direction, row vector
   k = 0L
   lambda = 0.001
   window,12,xsize=600,ysize=400, $
      title='FFN(scaled conjugate gradient)'
   wset,12
   progressbar = Obj_New('cgprogressbar', $
     Color='blue', Text='0', $
     title='Training: epoch No...',xsize=250,ysize=20)
   progressbar->start
   eivminmax = '?'
   repeat begin
      if progressbar->CheckCancel() then begin
         print,'Training interrupted'
         progressbar->Destroy
         return
      endif
      d2 = total(d*d)               ; d^2
      dTHd = total(self->Rop(d)*d)  ; d^T.H.d
      delta = dTHd+lambda*d2
      if delta lt 0 then begin
         lambda = 2*(lambda-delta/d2)
         delta = -dTHd
      endif
      E1 = self->cost(key)          ; E(w)
      (*self.cost_array)[k] = E1
      if key then $
          (*self.valid_cost_array)[k]=self->cost(2)
      dTg = total(d*g)              ; d^T.g
      alpha = -dTg/delta
      dw = alpha*d
      w = w+dw
      *self.Wh = reform(w[0:self.LL*(self.NN+1)-1], $
         self.LL,self.NN+1)
      *self.Wo = reform(w[self.LL*(self.NN+1):*], $
         self.KK,self.LL+1)
      E2 = self->cost(key)          ; E(w+dw)
      Ddelta = -2*(E1-E2)/(alpha*dTg); quadricity
      if Ddelta lt 0.25 then begin
         w = w - dw   ; undo weight change
         *self.Wh = reform(w[0:self.LL*(self.NN+1)-1],$
            self.LL,self.NN+1)
         *self.Wo = reform(w[self.LL*(self.NN+1):*], $
            self.KK,self.LL+1)
         lambda = 4*lambda       ; decrease step size
         if lambda gt 1e20 then $; if step too small
           k=self.epochs   $     ;  then give up
         else d = -g             ;  else restart
      end else begin
         k++
         if Ddelta gt 0.75 then lambda = lambda/2
         g = self->gradient()
         if k mod nw eq 0 then begin
             beta = 0
             eivs = self->eigenvalues()
             eivminmax = string(min(eivs)/max(eivs),$
                format='(F10.6)')
         end else beta = total(self->Rop(g)*d)/dTHd
         d = beta*d-g
         if k mod 10 eq 0 then plot,$
            *self.cost_array,xrange=[0,k>100L], $
            color=0, background='FFFFFF'XL, $
            ytitle='cross entropy', $
            xtitle='Epoch ['+ $
             'min(lambda)/max(lambda)='+ $
             eivminmax+']'
         if key then oplot, $
            *self.valid_cost_array,color=0,linestyle=2
      endelse
      progressbar->Update,k*100/self.epochs
   endrep until k ge self.epochs
   progressbar->Destroy
End

Pro FFNCG__Define
   class  = { FFNCG, $
              cost_array: ptr_new(), $
              valid_cost_array: ptr_new(), $
              GTs: ptr_new(), $ ; for training
                                ; (as opposed to validation)
              LTs: ptr_new(), $ ; ditto
              ntp: 0L, $        ; ditto
              epochs: 0L, $
              Inherits FFN $
            }
End