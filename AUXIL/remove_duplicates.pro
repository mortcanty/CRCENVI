; docformat = 'rst'
; remove_duplicates.pro
;+
; :Description:
;       remove duplicates in an integer array
; :Params:
;       array:  in, required, type = array of integer 
; :Returns:
;       array of integer                      
;-
Function remove_duplicates, array
   return, where( histogram(array,omin=om) gt 0)+om
End