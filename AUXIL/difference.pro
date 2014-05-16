; docformat = 'rst'
; difference.pro
;+
; :Description:
;       integer set difference
; :Params:
;       set1:
;       set2:  in, required, type = array of integer 
;          sets to be differenced
; :Returns:
;       array of integer 
; :Uses:
;       REMOVE_DUPLICATES                           
;-
Function difference, set1, set2
   return, where( histogram([remove_duplicates(set1),remove_duplicates(set2)],omin=om) eq 1 ) + om
End