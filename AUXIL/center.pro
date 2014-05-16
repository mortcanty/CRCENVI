; docformat = 'rst'
; center.pro
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
;       center a kernel matrix::
;          Shawe-Taylor, J. and Cristianini, N. (2004).
;          Kernel Methods for Pattern Analysis.
;          Cambridge University Press. 
; :Params:
;       K:  in, required, type = array of float 
;          kernel matrix to be centered
; :Returns:
;       array of float                 
; :Author:
;       Mort Canty (2009)       
;-
function center, K
   m = (size(K))[1]
   ones=dblarr(m,m)+1
   return, K-(ones##K+K##ones-total(K)/m)/m
end