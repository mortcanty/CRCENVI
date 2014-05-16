; docformat = 'rst'
; ffn3__define.pro
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
;     This is a generic class with no training method   
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
; :Author:      
;       Mort Canty (2009)            
Function FFN3::Init, Gs, Ls, L1, L2

common random_seed, sd

;sd=12345L

   catch, theError
   if theError ne 0 then begin
      catch, /cancel
      ok = dialog_message(!Error_State.Msg + ' Returning...', /error)
      return, 0
   endif
; network architecture
   self.LL1 = L1
   self.LL2 = L2
   self.np = n_elements(Gs[*,0])
   self.NN = n_elements(Gs[0,*])
   self.KK = n_elements(Ls[0,*])
; biased output vector from hidden layer (column vector)
   self.N1= ptr_new(dblarr(L1+1))
   self.N2= ptr_new(dblarr(L2+1))
; biased exemplars (column vectors)
   self.Gs = ptr_new([[dblarr(self.np)+1],[Gs]])
   self.Ls = ptr_new(Ls)
; weight matrices (each column is a neuron weight vector)
   self.W1h = ptr_new(randomu(sd,L1,self.NN+1,/double)-0.5)
   self.W2h = ptr_new(randomu(sd,L2,L1+1,/double)-0.5)
   self.Wo  = ptr_new(randomu(sd,self.KK,L2+1,/double)-0.5)
   return,1
End

Pro FFN3::Cleanup
   ptr_free, self.Gs
   ptr_free, self.Ls
   ptr_free, self.W1h
   ptr_free, self.W2h
   ptr_free, self.Wo
   ptr_free, self.N1
   ptr_free, self.N2
End

Function FFN3::forwardPass, x
; logistic activation for hidden neurons, N1,N2 set as side effect
   *self.N1 = [[1.0D],[1/(1+exp(-transpose(*self.W1h)##x))]]
   *self.N2 = [[1.0D],[1/(1+exp(-transpose(*self.W2h)##*self.N1))]]
; softmax activation for output neurons
   I = transpose(*self.Wo)##*self.N2
   A = exp(I-max(I))
   return, A/total(A)
End

Function FFN3::vForwardPass, Gs, N1, N2
; vectorized forward pass, returns output vector
; N1,N2 returned as a side effects
   np = n_elements(Gs[*,0])
   Ones = dblarr(np) + 1.0D
   N1 = [[Ones],[1/(1+exp(-transpose(*self.W1h)##Gs))]]
   N2 = [[Ones],[1/(1+exp(-transpose(*self.W2h)##N1))]]
   Io = transpose(*self.Wo)##N2
   maxIo = max(Io,dimension=2)
   for k=0,self.KK-1 do Io[*,k]=Io[*,k]-maxIo
   A = exp(Io)
   sum = total(A,2)
   M = dblarr(np,self.KK)
   for k=0,self.KK-1 do M[*,k] = A[*,k]/sum
   return, M
End

Function FFN3:: classify, X, Probs, labels=labels
; vectorized class membership probabilities (X not biased)
   nx = n_elements(X[*,0])
   Ones = dblarr(nx) + 1.0
   N1 = [[Ones],[1/(1+exp(-transpose(*self.W1h)##[[Ones],[X]]))]]
   N2 = [[Ones],[1/(1+exp(-transpose(*self.W2h)##N1))]]
   Io = transpose(*self.Wo)##N2
   maxIo = max(Io,dimension=2)
   for k=0,self.KK-1 do Io[*,k]=Io[*,k]-maxIo
   A = exp(Io)
   sum = total(A,2)
   Probs = dblarr(nx,self.KK)
   for k=0,self.KK-1 do Probs[*,k] = A[*,k]/sum
; vectorized class memberships
   void = max(probs,labels,dimension=2)
   return, byte(labels/nx+1)
End

Function FFN3:: cost, key
   case key of
      0: indices = lindgen(self.np)                                           ; all training data
      1: void = where( (lindgen(self.np) mod 2) eq 0, complement=indices )    ; odd-numbered training data
      2: indices = where( (lindgen(self.np) mod 2) eq 0 )                     ; even-numbered training data
   endcase
   np = n_elements(indices)
   Gs=(*self.Gs)[indices,*]
   Ls=(*self.Ls)[indices,*]
   Ones = dblarr(np) + 1.0D
   N1 = [[Ones],[1/(1+exp(-transpose(*self.W1h)##[Gs]))]]
   N2 = [[Ones],[1/(1+exp(-transpose(*self.W2h)##[N1]))]]
   Io = transpose(*self.Wo)##N2
   maxIo = max(Io,dimension=2)
   for k=0,self.KK-1 do Io[*,k]=Io[*,k]-maxIo
   A = exp(Io)
   sum = total(A,2)
   Ms = dblarr(np,self.KK)
; softmax
   for k=0,self.KK-1 do Ms[*,k] = A[*,k]/sum
   return, -total(Ls*alog(Ms+1e-20))
End

Pro FFN3__Define
class =  { FFN3, $
           NN: 0L,               $ ;input dimension
           LL1: 0L,              $ ;number of 1st layer hidden units
           LL2: 0L,              $ ;number of 2nd layer hidden units
           KK: 0L,               $ ;output dimension
           W1h:ptr_new(),        $ ;1st layer hidden weights
           W2h:ptr_new(),        $ ;2nd layer hidden weights
           Wo:ptr_new(),         $ ;output weights
           Gs:ptr_new(),         $ ;training pairs
           Ls:ptr_new(),         $
           N1:ptr_new(),         $ ;output vector from 1st hidden layer
           N2:ptr_new(),         $ ;output vector from 2nd hidden layer
           np: 0L                $ ;number of training pairs
          }
End