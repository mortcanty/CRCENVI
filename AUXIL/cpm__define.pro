; docformat = 'rst'
; cpm__define.pro
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
;       Object class for sequential variance-
;       covariance matrix calculation using 
;       the method of provisional means     
; :Params:             
;     p: in,required
;        dimension of observations       
; :Examples:         
;        cpm = Obj_New("CPM",p) 
;        
;        cpm->update(Xs,weights=Ws) 
; :Uses:
;     provmeans.dll
;     provmeans.dlm                           
; :Author:      
;       Mort Canty (2009)  
;-
Function CPM::Init, NN
   self.mn = ptr_new(dblarr(NN))
   self.cov= ptr_new(dblarr(NN,NN))
   self.sw = 0.0000001D
   return, 1
End

Pro CPM::Cleanup
   Ptr_Free, self.mn
   Ptr_Free, self.cov
End

Pro CPM::Update, Xs, weights = Ws
   COMPILE_OPT STRICTARR
   sz = size(Xs)
   NN = sz[1]
   n = sz[2]
   if n_elements(Ws) eq 0 then Ws = fltarr(sz[2])+1.0
   sw  =  self.sw
   mn  = *self.mn
   cov = *self.cov
   void = provmeans(float(Xs),float(Ws),NN,n,sw,mn,cov)          
   self.sw   = sw
   *self.mn  = mn
   *self.cov = cov
End

Function CPM::Covariance
   c = *self.cov/(self.sw-1)
   d = diag_matrix(diag_matrix(c))
   return, c+transpose(c)-d
End

Function CPM::Means
   return, *self.mn
End

Pro CPM__Define
class = {CPM, $
         sw: 0.0D, $     ;current sum of weights
         mn: ptr_new(), $;current mean
         cov: ptr_new()} ;current sum of cross products
End