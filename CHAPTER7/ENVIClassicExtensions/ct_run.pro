; docformat = 'rst'
; ct_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO ct_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Contingency Table', $
      REF_VALUE = 'Overlay Classes', $
      EVENT_PRO = 'ct_run', $
      UVALUE = 'CT',$
      POSITION = 'after'
END

;+
; :Description:
;       ENVI extension to determine contingency
;       table (confusion matrix) and classification
;       accuracies from test classifcation results::
;          Richards, J. A. and Jia, X. (2006). Remote 
;          Sensing Digital Image Analysis ;(4th Ed.). 
;          Springer.       
; :Params:
;       event:  in, required 
;          if called from ENVI                  
; :Uses:
;       ENVI
; :Author:
;       Mort Canty (2013)     
;-
Pro CT_Run, event

COMPILE_OPT IDL2

   print, '---------------------------------'
   print, 'Classification Statistics'
   print, systime(0)
   print, '---------------------------------'

   fileName = Dialog_Pickfile(Filter=['*.tst'],/Read)
   if fileName eq '' then begin
      print, 'cancelled'
      return
   endif

   catch, theError
   if theError ne 0 then begin
      print, 'Error reading file: '+fileName
      return
   endif
   
   openR, lun, fileName, /get_lun
; read and print the 4 line header
   aLine = ''
   for i=0,3 do begin
      readF, lun, aLine
      print, aLine
   endfor
; get number of test data and classes
   n = 0L & K = 0
   readF, lun, n, K
   print, 'Test observations: '+strtrim(n,2)
   print, 'Classes: '+strtrim(K,2)
; contingency table
   CT = fltarr(K+2,K+2)
; fill the CT:
   y = 0L
   k1 = 0
   k2 = 0
   y = 0L
   for i=0,n-1 do begin
      readF, lun, k1, k2
      (CT[k2-1,k1-1])++
      if k1 ne k2 then y++
   endfor
   free_lun, lun
; fill the contingency table
   CT[*,K] = total(CT,2)
   CT[K,*] = total(CT,1)
   for i=0,K-1 do begin
      CT[i,K+1] = CT[i,i]/CT[i,K]
      CT[K+1,i] = CT[i,i]/CT[K,i]
   endfor
; overall misclassification rate
   sigma = sqrt((float(y)*(n-y))/(float(n)*n*n))
   low = (y+1.921-1.96*sqrt(0.96+y*(n-y)/float(n)))/(3.842+n)
   high= (y+1.921+1.96*sqrt(0.96+y*(n-y)/float(n)))/(3.842+n)
   print, 'Misclassification rate',float(y)/n, format='(A23,F8.4)'
   print, 'Standard deviation', sigma, format='(A23,F8.4)'
   print, 'Conf. interval (95%)', low, high, format='(A23,2F8.4)'
; Kappa coefficient
   t1 = float(n-y)/n
   t2 = total(CT[0:K-1,K]*transpose(CT[K,0:K-1]))/(float(n)*n)
   Kappa = (t1 - t2)/(1 - t2)
   t3 = 0.0
   for i=0,K-1 do t3 = t3 + CT[i,i]*(CT[i,K]+CT[K,i])
   t3 = t3/(float(n)*n)
   t4 = 0.0
   for i=0,K-1 do for j=0,K-1 do t4 = t4 + CT[i,j]*(CT[j,K]+CT[K,i])^2
   t4 = t4/(float(n)*n*n)
   sigma2 = t1*(1-t1)/(1-t2)^2
   sigma2 = sigma2 + 2*(1-t1)*(2*t1*t2-t3)/(1-t2)^3
   sigma2 = sigma2 + ((1-t1)^2)*(t4-4*t2^2)/(1-t2)^4
   sigma = sqrt(sigma2/n)
   print, 'Kappa coefficient',Kappa, format='(A23,F8.4)'
   print, 'Standard deviation', sigma, format='(A23,F8.4)'
   print, 'Contingency Table'
   format = '('+strtrim(K+1,2)+'I7,F7.3)'
   for i=0,K do print, CT[*,i],format=format
   format = '('+strtrim(K+2,2)+'F7.3)'
   print, CT[*,K+1], format=format

End