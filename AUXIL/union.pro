; docformat = 'rst'
; union.pro
;+
; :Description:
;       integer set union
; :Params:
;       set1:
;       set2:  in, required, type = array of integer 
;          sets to be combined
; :Returns:
;       array of integer 
; :Uses:
;       REMOVE_DUPLICATES                           
;-
Function union, set1, set2
   return, remove_duplicates([set1,set2])
End