; docformat = 'rst'
; complement.pro
;+
; :Description:
;       integer set complement
; :Params:
;       set:  in, required, type = array of integer 
;          set to be complemented
; :Returns:
;       array of integer 
; :Uses:
;       DIFFERENCE                           
;-
Function complement, set
   return, difference(set,lindgen(max(set)+1))
End