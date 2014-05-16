; docformat = 'rst'
; ffnkal__define.pro
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
;     Implements Kalman filter training::
;        Shah, S. and Palmieri, F. (1990). Meka — A fast, 
;        local algorithm for training feed forward neural 
;        networks. Proceedings of the International Joint
;        Conference on Neural Networks, San Diego, I(3), 
;        41–46.
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
;     ffn = Obj_New("FFNKAL",Gs,Ls,L)
; :Uses:
;       COYOTE 
; :Author:      
;       Mort Canty (2009)            
;-
Function FFNKAL::Init, Gs, Ls, L
   catch, theError
   if theError ne 0 then begin
      catch, /cancel
      ok = dialog_message(!Error_State.Msg + ' Returning...', /error)
      return, 0
   endif
; initialize the superclass
   if not self->FFN::Init(Gs, Ls, L) then return, 0
   self.iterations = 1000*self.np
   self.cost_array = ptr_new(dblarr((self.iterations+100)/100))
   self.valid_cost_array = ptr_new(dblarr((self.iterations+100)/100))
; weight covariance matrices for hidden neurons
   self.Sh = ptr_new(dblarr(self.NN+1,self.NN+1,L))
   for i=0,L-1 do (*self.Sh)[*,*,i]=diag_matrix(dblarr(self.NN+1)+100)
; weight covariance matrices for output neurons
   self.So = ptr_new(dblarr(L+1,L+1,self.KK))
   for i=0,self.KK-1 do (*self.So)[*,*,i]=diag_matrix(dblarr(L+1)+100)
   return, 1
End

Pro FFNKAL::Cleanup
   ptr_free, self.cost_array
   ptr_free, self.valid_cost_array
   ptr_free, self.Sh
   ptr_free, self.So
   self->FFN::Cleanup
End

Pro FFNKAL:: train, key=key
   if n_elements(key) eq 0 then key=0
; define update matrices for Wh and Wo
   dWh = dblarr(self.LL,self.NN+1)
   dWo = dblarr(self.KK,self.LL+1)
   iter = 0L
   iter100 = 0L
   progressbar = Obj_New('cgprogressbar', $
     Color='blue', Text='0',$
     title='Training: example number...', $
     xsize=250,ysize=20)
   progressbar->start
   window,12,xsize=600,ysize=400, $
     title='FFN(Kalman filter)'
   wset,12
   repeat begin
      if progressbar->CheckCancel() then begin
         print,'Training interrupted'
         progressbar->Destroy
         return
      endif
; select training pair at random
      ell = long(self.np*randomu(seed))
; use only odd ones if validation option is set
      if key eq 1 then ell = (ell*2 mod (self.np-2))+1
      x=(*self.Gs)[ell,*]
      y=(*self.Ls)[ell,*]
; send it through the network
      m=self->forwardPass(x)
; error at output
      e=y-m
; loop over the output neurons
      for k=0,self.KK-1 do begin
;    linearized input (column vector)
         Ao = m[k]*(1-m[k])*(*self.N)
;    Kalman gain
         So = (*self.So)[*,*,k]
         SA = So##Ao
         Ko = SA/((transpose(Ao)##SA)[0]+1)
;    determine delta for this neuron
         dWo[k,*] = Ko*e[k]
;    update its covariance matrix
         So = So - Ko##transpose(Ao)##So
         (*self.So)[*,*,k] = So
      endfor
; update the output weights
      *self.Wo = *self.Wo + dWo
; backpropagated error
      beta_o =e*m*(1-m)
; loop over the hidden neurons
      for j=0,self.LL-1 do begin
;   linearized input (column vector)
         Ah = X*(*self.N)[j+1]*(1-(*self.N)[j+1])
;   Kalman gain
         Sh = (*self.Sh)[*,*,j]
         SA = Sh##Ah
         Kh = SA/((transpose(Ah)##SA)[0]+1)
;   determine delta for this neuron
         dWh[j,*] = Kh*((*self.Wo)[*,j+1]##beta_o)[0]
;   update its covariance matrix
         Sh = Sh - Kh##transpose(Ah)##Sh
         (*self.Sh)[*,*,j] = Sh
   endfor
; update the hidden weights
      *self.Wh = *self.Wh + dWh
; record cost history
      if iter mod 100 eq 0 then begin
         if key then begin
            (*self.cost_array)[iter100]= $
             alog10(self->cost(1)) ; training cost
            (*self.valid_cost_array)[iter100]= $
             alog10(self->cost(2)) ; validation cost
         end else (*self.cost_array)[iter100]= $
             alog10(self->cost(0)) ; full training cost
         iter100 = iter100+1
         progressbar->Update,iter*100/self.iterations
         plot,*self.cost_array,xrange=[0,iter100], $
            color=0,background='FFFFFF'XL,$
            xtitle='Iterations/100)', $
            ytitle='log(cross entropy)'
         if key then oplot, *self.valid_cost_array, $
            color=0,linestyle=2
      end
      iter=iter+1
   endrep until iter eq self.iterations
   progressbar->Destroy
End

Pro FFNKAL__Define
   class  = { FFNKAL, $
              cost_array: ptr_new(), $
              valid_cost_array: ptr_new(), $
              iterations: 0L, $
              Sh: Ptr_New(), $
              So: Ptr_New(), $
              Inherits FFN $
            }
End