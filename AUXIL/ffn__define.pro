; docformat = 'rst'
; ffn__define.pro
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
;     This is a generic class with no training method        
; :Params:
;        Gs: in, required 
;           array of observation column vectors     
;        Ls: in, required
;            array of class label column vectors
;            of form (0,0,1,0,0,...0)^T          
;        L: in, required
;           number of hidden neurons
; :Author:      
;       Mort Canty (2009)           
;-
Function FFN::Init, Gs, Ls, L
; network architecture
   self.LL = L
   self.np = n_elements(Gs[*,0])
   self.NN = n_elements(Gs[0,*])
   self.KK = n_elements(Ls[0,*])
; biased output vector from hidden layer 
   self.N= ptr_new(dblarr(L+1))
; biased exemplars (column vectors)
   self.Gs = ptr_new([[dblarr(self.np)+1],[Gs]])
   self.Ls = ptr_new(Ls)
; weight matrices 
   self.Wh = ptr_new(randomu(seed,L,self.NN+1,/double)$
                                                -0.5)
   self.Wo = ptr_new(randomu(seed,self.KK,L+1,/double)$
                                                -0.5)
   return,1
End

Function FFN::forwardPass, G
; forward pass through network 
   expnt = transpose(*self.Wh)##G
   *self.N = [[1.0],[1/(1+exp(-expnt))]]
; softmax activation for output neurons
   I = transpose(*self.Wo)##*self.N
   A = exp(I-max(I))
   return, A/total(A)
End

Function FFN::classify, Gs, Probs
; vectorized class membership probabilities
   nx = n_elements(Gs[*,0])
   Ones = dblarr(nx) + 1.0
   expnt = transpose(*self.Wh)##[[Ones],[Gs]]
   N = [[Ones],[1/(1+exp(-expnt))]]
   Io = transpose(*self.Wo)##N
   maxIo = max(Io,dimension=2)
   for k=0,self.KK-1 do Io[*,k]=Io[*,k]-maxIo
   A = exp(Io)
   sum = total(A,2)
   Probs = fltarr(nx,self.KK)
   for k=0,self.KK-1 do Probs[*,k] = A[*,k]/sum
   void = max(probs,labels,dimension=2)
   return, byte(labels/nx+1)
End

Function FFN::vForwardPass, Gs, N
; vectorized forward pass
   np = n_elements(Gs[*,0])
   Ones = dblarr(np) + 1.0
   expnt = transpose(*self.Wh)##Gs
   N = [[Ones],[1/(1+exp(-expnt))]]
   Io = transpose(*self.Wo)##N
   maxIo = max(Io,dimension=2)
   for k=0,self.KK-1 do Io[*,k]=Io[*,k]-maxIo
   A = exp(Io)
   sum = total(A,2)
   M = dblarr(np,self.KK)
   for k=0,self.KK-1 do M[*,k] = A[*,k]/sum
   return, M
End

Function FFN:: cost, key
   case key of
      0: indices = lindgen(self.np)                                           ; all training data
      1: void = where( (lindgen(self.np) mod 2) eq 0, complement=indices )    ; odd-numbered training data
      2: indices = where( (lindgen(self.np) mod 2) eq 0 )                     ; even-numbered training data
   endcase
   np = n_elements(indices)
   Gs=(*self.Gs)[indices,*]
   Ls=(*self.Ls)[indices,*]
   Ones = dblarr(np) + 1.0
   expnt = transpose(*self.Wh)##[Gs]
   N = [[Ones],[1/(1+exp(-expnt))]]
   Io = transpose(*self.Wo)##N
   maxIo = max(Io,dimension=2)
   for k=0,self.KK-1 do Io[*,k]=Io[*,k]-maxIo
   A = exp(Io)
   sum = total(A,2)
   Ms = dblarr(np,self.KK)
; softmax
   for k=0,self.KK-1 do Ms[*,k] = A[*,k]/sum
   return, -total(Ls*alog(Ms+1e-20))
End

Function FFN::GetWeights
   return, [(*self.Wh)[*],(*self.Wo)[*]]
End

Pro FFN::PutWeights, w
   *self.Wh = reform(w[0:self.LL*(self.NN+1)-1],self.LL,self.NN+1)
   *self.Wo = reform(w[self.LL*(self.NN+1):*],self.KK,self.LL+1)
End

Pro FFN::Cleanup
   ptr_free, self.Gs
   ptr_free, self.Ls
   ptr_free, self.Wh
   ptr_free, self.Wo
   ptr_free, self.N
End

Pro FFN__Define
class =  { FFN, $
           NN: 0L,               $ ;input dimension
           LL: 0L,               $ ;number of hidden units
           KK: 0L,               $ ;output dimension
           Wh:ptr_new(),         $ ;hidden weights
           Wo:ptr_new(),         $ ;output weights
           Gs:ptr_new(),         $ ;training pairs
           Ls:ptr_new(),         $
           N:ptr_new(),          $ ;output vector from hidden layer
           np: 0L                $ ;number of training pairs
          }
End