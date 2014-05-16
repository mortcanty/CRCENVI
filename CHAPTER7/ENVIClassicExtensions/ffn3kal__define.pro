; docformat = 'rst'
; ffn3kal__define.pro
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
;       Object class for implementation of a
;       three-layer, feed-forward neural network
;       for classification of multi-spectral images
;       with AdaBoost.M1 or .M2 boosting
;       
;       Implements Kalman filter training.
;       
;       Refs: Shaw and Palmieri, Proc. Int. Joint Conf.
;           on Neural Networks, San Diego (1990) I(3), 41-46
;           
;           Freund and Shapire, J. Comp. and Syst. Sci. 55, 1997, 119-139
;           Schwenk and Bengio, Neural Computation 12(8) 2000, 1869-1887
;    
;  ffn = Obj_New("FFN3KAL",Xs,Ls,L1,L2,D=D,epochs=epochs,M2=M2)
; 
;       Xs:      array of observation column vectors
;                observations are in [0,1]
;                
;       Ls:      array of class label column vectors
;                of form (0,0,1,0,0,...0)^T
;                
;       L1,L2:   numbers of hidden neurons
; 
;       D:       adaboost sample weight distribution
;       
;       epochs:  numberof training epochs
;       
;       M2:      if 0, then AdaBoost.M1 (default), else M2
; :Uses:
;       PROGRESSBAR 
; :Author:      
;       Mort Canty (2009)
;       Juelich Research Center
;       m.canty@fz-juelich.de              
;-

; docformat = 'rst'
; ffn3kal__define.pro
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
;     three-layer, feed-forward neural network classifier.
;     Implements Kalman filter training::
;        Shah, S. and Palmieri, F. (1990). Meka — A fast, 
;        local algorithm for training feed forward neural 
;        networks. Proceedings of the International Joint
;        Conference on Neural Networks, San Diego, I(3), 
;        41–46.
; :Inherits:
;        FFN3__DEFINE       
; :Params:
;        Gs: in, required 
;           array of observation column vectors     
;        Ls: in, required
;            array of class label column vectors
;            of form (0,0,1,0,0,...0)^T          
;        L1, L2: in, required
;           number of hidden neurons
; :Keywords:
;        epochs: in, optional
;           training epochs per network in emsemble
;        D: in, required 
;           training sample distribution              
; :Examples:   
;     ffn = Obj_New("FFN3KAL",Gs,Ls,L1,L2)
; :Uses:
;       COYOTE 
; :Author:      
;       Mort Canty (2009)            
Function FFN3KAL::Init, Gs, Ls, L1, L2, epochs=epochs, D=D
   catch, theError
   if theError ne 0 then begin
      catch, /cancel
      ok = dialog_message(!Error_State.Msg + ' Returning...', /error)
      return, 0
   endif
; initialize the superclass
   if not self->FFN3::Init(Gs, Ls, L1, L2) then return, 0
   if n_elements(D) eq 0 then D = dblarr(self.np)+1.0/self.np
   self.D = ptr_new(D)
; cumulative sample distribution for random selection
   tmp = shift(total(*self.D,/cumulative),1)
   tmp[0]=0.0
   self.SD = ptr_new(tmp)
; interations
   if n_elements(epochs) eq 0 then epochs = 1
   self.iterations = epochs*self.np
   self.cost_array = ptr_new(dblarr((self.iterations+100)/100))
; weight covariance matrices for hidden neurons
   self.S1h = ptr_new(dblarr(self.NN+1,self.NN+1,L1))
   self.S2h = ptr_new(dblarr(L1+1,L1+1,L2))
   for i=0,L1-1 do (*self.S1h)[*,*,i]=diag_matrix(dblarr(self.NN+1)+100)
   for i=0,L2-1 do (*self.S2h)[*,*,i]=diag_matrix(dblarr(L1+1)+100)
; weight covariance matrices for output neurons
   self.So = ptr_new(dblarr(L2+1,L2+1,self.KK))
   for i=0,self.KK-1 do (*self.So)[*,*,i]=diag_matrix(dblarr(L2+1)+100)
   return, 1
End

Pro FFN3KAL::Cleanup
   ptr_free, self.cost_array
   ptr_free, self.S1h
   ptr_free, self.S2h
   ptr_free, self.So
   ptr_free, self.SD
   ptr_free, self.D
   ptr_free, self.V
   self->FFN3::Cleanup
End


Pro FFN3KAL::train, verbose=verbose, cancel=cancel
   on_error, 2
   cancel = 0
   if n_elements(verbose) eq 0 then verbose = 0
; update matrices for W1h, W2h and Wo
   dW1h = dblarr(self.LL1,self.NN+1)
   dW2h = dblarr(self.LL2,self.LL1+1)
   dWo =  dblarr(self.KK, self.LL2+1)
   work = dblarr(self.LL2)
   iter = (iter100 = 0L)
   progressbar = Obj_New('progressbar', $
     Color='blue', Text='0',$
     title='Training: example number...', $
     xsize=250,ysize=20)
   progressbar->start
   if verbose then begin
      window,12,xsize=600,ysize=400, $
        title='FFN(3-layer, Kalman filter)'
      wset,12
   endif
   repeat begin
      if progressbar->CheckCancel() then begin
         print,'Training interrupted'
         cancel = 1
         progressbar->Destroy
         return
      endif
; select exemplar pair at random from sample distribution
      ell = long(total(*self.SD lt randomu(seed))-1)
      x=(*self.Gs)[ell,*]
      y=(*self.Ls)[ell,*]
; send it through the network
      m=self->forwardPass(x)
; error at output
      e=(y-m)
      mm = m*(1-m)
      nn2= (*self.N2)*(1-(*self.N2))
      nn1= (*self.N1)*(1-(*self.N1))
; loop over the output neurons
      for k=0,self.KK-1 do begin
;    linearized input (column vector)
         Ao = (*self.N2)*mm[k]
;    Kalman gain
         So = (*self.So)[*,*,k]
         SA = So##Ao
         Ko = SA/((transpose(Ao)##SA)[0]+1)
;    determine delta for this neuron
         dWo[k,*] = Ko*e[k]
;    update its covariance matrix
         (*self.So)[*,*,k] = (So -= Ko##transpose(Ao)##So)
      endfor
; update the output weights
      *self.Wo += dWo
; backpropagated error
      beta_o =e*mm
; loop over the hidden neurons in 2nd layer
      for j=0,self.LL2-1 do begin
;   linearized input (column vector)
         Ah = (*self.N1)*nn2[j+1]
;   Kalman gain
         Sh = (*self.S2h)[*,*,j]
         SA = Sh##Ah
         Kh = SA/((transpose(Ah)##SA)[0]+1)
;   determine delta for this neuron
         work[j] = ((*self.Wo)[*,j+1]##beta_o)[0]
         dW2h[j,*] = Kh*work[j]
;   update its covariance matrix
         (*self.S2h)[*,*,j] = (Sh -= Kh##transpose(Ah)##Sh)
      endfor
; update the 2nd layer hidden weights
      *self.W2h += dW2h
; loop over the hidden neurons in 1st layer
      for j=0,self.LL1-1 do begin
;   linearized input (column vector)
         Ah = X*nn1[j+1]
;   Kalman gain
         Sh = (*self.S1h)[*,*,j]
         SA = Sh##Ah
         Kh = SA/((transpose(Ah)##SA)[0]+1)
;   determine delta for this neuron
         N2 = (*self.N2)[1:self.LL2]
         dW1h[j,*] = Kh*total( work*((*self.W2h)[*,j+1])*N2*(1-N2) )
;   update its covariance matrix
         (*self.S1h)[*,*,j] = (Sh -= Kh##transpose(Ah)##Sh)
      endfor
; update the 1st layer hidden weights
      *self.W1h += dW1h
; record cost history
      if iter mod 100 eq 0 then begin
         if verbose  then (*self.cost_array)[iter100]= $
             self->cost() ; full training cost
         iter100 = iter100+1
         progressbar->Update,iter*100/self.iterations,$
               text=strtrim(iter,2)
         if verbose then plot,*self.cost_array,xrange=[0,iter100], $
            color=0,background='FFFFFF'XL,$
            xtitle='Iterations/100)', $
            ytitle='Quadratic error'
      end
      iter=iter+1
   endrep until iter eq self.iterations
   progressbar->Destroy
End

Function FFN3KAL:: cost
   Gs=(*self.Gs)
   Ls=(*self.Ls)
   Ones = dblarr(self.np) + 1.0D
   N1 = [[Ones],[1/(1+exp(-transpose(*self.W1h)##[Gs]))]]
   N2 = [[Ones],[1/(1+exp(-transpose(*self.W2h)##[N1]))]]
   Io = transpose(*self.Wo)##N2
   maxIo = max(Io,dimension=2)
   for k=0,self.KK-1 do Io[*,k]=Io[*,k]-maxIo
   A = exp(Io)
   sum = total(A,2)
   Ms = dblarr(self.np,self.KK)
; softmax
   for k=0,self.KK-1 do Ms[*,k] = A[*,k]/sum
   return, 0.5*total(total((Ls-Ms)^2,2)*(*self.D))
End

Pro FFN3KAL__Define
   class  = { FFN3KAL, $
              cost_array: ptr_new(), $
              SD:ptr_new(),   $ ;training sample distribution;
              D: ptr_new(),   $ ;sample weights
              V: ptr_new(),   $ ;cost function weights
              iterations: 0L, $
              S1h: Ptr_New(), $ ;weight covariance matrices
              S2h: Ptr_New(), $
              So:  Ptr_New(), $
              Inherits FFN3   }
End