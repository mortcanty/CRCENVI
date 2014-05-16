; docformat = 'rst'
; ortho_regress.pro
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
;      Orthogonal regression between two vectors::    
;        Canty, M. J., Nielsen, A. A., and Schmidt, M. 
;        (2004). Automatic radiometric normalization 
;        of multitemporal satellite imagery. 
;        Remote Sensing of Environment, 
;        91(3-4), 441–451. I
; :Params:
;       X: in,required
;       Y: in,required
;       b: out, required 
;       Xm: out, required 
;       Ym: out, required 
;       sigma_a: out, required        
;       sigma_b: out, required        
;       sigma(RMSE): out, required   
; :Keywords:
;       rank: out,optional
;           Spearman rank order correlation coefficient
;           (two-element vector containing the rank 
;           correlation coefficient and the two-sided 
;           significance of its deviation from zero)           
; :Author:
;       Mort Canty (2009)      
;-
pro ortho_regress, X, Y, b, Xm, Ym, $
                   sigma_a, sigma_b, sigma, rank=rank
   m = n_elements(X)
   Xm = mean(X)
   Ym = mean(Y)
   S = correlate([X,Y],/covariance,/double)
   Sxx = S[0,0]
   Syy = S[1,1]
   Sxy = S[1,0]
   void = eigenql(S,eigenvectors=eigenvectors,/double)
; slope
   b = eigenvectors[1,0]/eigenvectors[0,0]
; standard errors
   sigma2 = m*(Syy-2*b*Sxy+b*b*Sxx)/((m-2)*(1+b^2))
   tau = sigma2*b/((1+b^2)*Sxy)
   sigma_b=sqrt( sigma2*b*(1+b^2)*(1+tau)/(m*Sxy))
   sigma_a=sqrt( sigma2*b*(1+b^2)*(Xm^2*(1+tau)+Sxy/b)$
                                          /(m*Sxy) )
   sigma = sqrt(sigma2)
   rank = r_correlate(X,Y)
end