; docformat = 'rst'
; intersection.pro
;+
; :Description:
;       integer set intersection
; :Params:
;       set1:
;       set2:  in, required, type = array of integer 
;          sets to be intersected
; :Returns:
;       array of integer                   
;-
Function intersection, set1, set2
   res =  where( histogram(set1,omin=om) gt 0 and histogram(set2,min=om) gt 0, count) + om
   if count gt 0 then return, res else return, -1
End