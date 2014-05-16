; docformat = 'rst'
; mcnemar_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO mcnemar_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'McNemar Test', $
      REF_VALUE = 'Overlay Classes', $
      EVENT_PRO = 'mcNemar_run', $
      UVALUE = 'MCNEMAR',$
      POSITION = 'after'
END

;+
; :Description:
;       ENVI extension to compare classifiers
;       on the basis of misclassifications
;       using McNemar's statistic::
;          Dietterich, T. G. (1998). Approximate 
;          statistical  tests for comparing supervised 
;          classification learning algorithms. 
;          Neural Computation, 10(7), 1895–;1923.       
; :Params:
;       event:  in, required 
;          if called from ENVI                  
; :Uses:
;       ENVI
; :Author:
;       Mort Canty (2013)  
;-
Pro McNemar_Run, event

COMPILE_OPT STRICTARR

   print, '--------------------------'
   print, 'Classification Comparison'
   print, systime(0)
   print, '--------------------------'

   fileNameA = Dialog_Pickfile(Filter=['*.tst'],/Read,title='Select first test result file')
   if fileNameA eq '' then begin
      print, 'cancelled'
      return
   endif

   fileNameB = Dialog_Pickfile(Filter=['*.tst'],/Read,title='Select second test result file')
   if fileNameB eq '' then begin
      print, 'cancelled'
      return
   endif

   !error_state.code=0
   catch, theError
   if theError ne 0 then begin
      print, 'Error reading file: '+fileName
      return
   endif
   openR, lunA, fileNameA, /get_lun
   openR, lunB, fileNameB, /get_lun
; read and print the 4 line headers
   print, 'First classifier'
   aLine = ''
   readF, lunA, aLine &  print, aLine
   readF, lunA, aLine &  print, aLine
   readF, lunA, aLine &  print, aLine
   readF, lunA, aLine &  print, aLine
   print, 'Second classifier'
   readF, lunB, aLine &  print, aLine
   readF, lunB, aLine &  print, aLine
   readF, lunB, aLine &  print, aLine
   readF, lunB, aLine &  print, aLine
; get number of test data and classes
   n = 0L & K = 0
   readF, lunA, n, K
   readF, lunB, nb, Kb
   if (n ne nb) or (K ne Kb) then begin
      print, 'Error: test files do not correspond'
      return
   endif
   print, 'Test observations: '+strtrim(n,2)
   print, 'Classes: '+strtrim(K,2)
; calculate McNemar
   y10 = 0.0
   y01 = 0.0
   for i=0,n-1 do begin
      readF, lunA, k1A, k2A
      readF, lunB, k1B, k2B
         if (k1A ne k2A) and (k1B eq k2B) then y10=y10+1
         if (k1A eq k2A) and (k1B ne k2B) then y01=y01+1
   endfor
   free_lun, lunA
   free_lun, lunB
   McN = (abs(y01-y10)-1.0)^2/(y01+y10)
   print, 'First classifier:', y10, format = '(A22,I8)'
   print, 'Second classifier:', y01, format = '(A22,I8)'
   print, 'McNemar statistic:', McN, format = '(A22,F8.4)'
   print, 'P-value:', 1.0-chisqr_pdf(McN,1), format = '(A22,F8.4)'
End