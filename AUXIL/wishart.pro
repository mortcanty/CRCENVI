; docformat = 'rst'
; wishart.pro
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
;       Test two multivariate distributions for
;       equal means and covariance matrices
;       assuming both parameter sets are estimated
;       with the same sample size.
;       Returns p-value::
;          Anderson, T. W. (2003). An Introduction 
;          to Multivariate Statistical Analysis.
;          Wiley Series in Probability and Statistics, 
;          third edition.
; :Params:
;       M1: in,required
;          mean of first distribution
;       M2: in,required
;          mean of second distribution
;       S1: in,required
;          variance of first distribution  
;       S2: in,required
;          variance of second distribution  
;       n: in,required
;          common number of samples
; :Keywords:         
;       statistic: out,optional
;          test statistic -2*rho*log(lambda)   
; :Author:
;       Mort Canty (2009)     
;-
FUNCTION Wishart, M1, M2, S1, S2, n, statistic=statistic
   NN = n_elements(M1)
   B = n*(S1 + S2 + 0.5*transpose(M1-M2)##(M1-M2))
   logLambda = n*(NN*alog(2*n)+0.5*alog(determ(S1))+0.5*alog(determ(S2))-alog(determ(B)))
   rho = 1.0-1.5*(2.0*NN^2+9.0*NN+11.0)/(n*(6.0*(NN+3.0)))
   w2  = NN*(NN+3.0)*1.5*((NN+1)*(NN+2)-6.0*(1.0-rho^2))/(n^2*48.0*rho^2)
   df  = NN*(NN+1)/2
   statistic = -2.0*rho*logLambda
   Pr1 = chisqr_pdf(statistic,df)
   Pr2 = chisqr_pdf(statistic,df+4)
   return, 1-(Pr1+w2*(Pr2-Pr1))
END