; docformat = 'rst'
; ffnbp__define.pro
; :Copyright:
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
;     Object class for implementation of a
;     two-layer, feed-forward neural network classifier.
;     Implements ordinary backpropagation training.
;    
;     ffn = Obj_New("FFNBP",Gs,Ls,L)
;
;     Gs:  array of observation column vectors
;     
;     Ls:  array of class label column vectors
;          of form (0,0,1,0,0,...0)^T
;          
;     L:   number of hidden neurons
; :Uses:
;       COYOTE 
; :Author:      
;       Mort Canty (2009)
;       Juelich Research Center
;       m.canty@fz-juelich.de              
;-

Function FFNBP::Init, Gs, Ls, L
   catch, theError
   if theError ne 0 then begin
      catch, /cancel
      ok = dialog_message(!Error_State.Msg + ' Returning...', /error)
      return, 0
   endif
; initialize the superclass
   if not self->FFN::Init(Gs, Ls, L) then return, 0
   self.cost_array = ptr_new(fltarr((100*self.np+100)/100))
   return, 1
End

Pro FFNBP::Cleanup
   ptr_free, self.cost_array
   self->FFN::Cleanup
End

Pro FFNBP::Train
   iter = 0L & epoch = 0L
   iterations = 100*self.np
   eta = 0.01 &  alpha = 0.5 ; learn rate & momentum
   progressbar=Obj_New('cgprogressbar', $
    /cancel,Text='0',$
    title='Training: example No...',xsize=250,ysize=20)
   progressbar->start
   window,12,xsize=600,ysize=400,title='Cost Function'
   wset,12
   inc_o1 = 0 & inc_h1 = 0
   repeat begin
      if progressbar->CheckCancel() then begin
         print,'Training interrupted'
         progressbar->Destroy & return
      endif
; select example pair at random
      nu = long(self.np*randomu(seed))
      x=(*self.Gs)[nu,*]
      ell=(*self.Ls)[nu,*]
; send it through the network
      m=self->forwardPass(x)
; determine the deltas
      d_o = ell - m
      d_h = (*self.N*(1-*self.N)* $
            (*self.Wo##d_o))[1:self.LL]
; update the synaptic weights
      inc_o = eta*(*self.N##transpose(d_o))
      inc_h = eta*(x##d_h)
      *self.Wo = *self.Wo + inc_o + alpha*inc_o1
      *self.Wh = *self.Wh + inc_h + alpha*inc_h1
      inc_o1 = inc_o & inc_h1 = inc_h
; record cost history
      if iter mod self.np eq 0 then begin
         (*self.cost_array)[epoch] = $
                    alog10(self->cost(0))
         epoch++
         progressbar->Update,iter*100/iterations
         plot,*self.cost_array,xrange=[0,epoch],$
             color=0,background='FFFFFF'XL,$
            xtitle='Epoch)',$
            ytitle='log(cross entropy)'
      end
      iter=iter+1
   endrep until iter eq iterations
   progressbar->destroy
End

Pro FFNBP__Define
   class  = { FFNBP, $
              cost_array: ptr_new(), $
              Inherits FFN $
            }
End