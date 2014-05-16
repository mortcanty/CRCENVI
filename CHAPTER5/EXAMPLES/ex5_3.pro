pro Ex5_3

; Airplane
   A = float([[0,0,0,0,0,1,0,0,0,0,0], $
             [0,0,0,0,1,1,1,0,0,0,0], $
             [0,0,0,0,1,1,1,0,0,0,0], $
             [0,0,0,1,1,1,1,1,0,0,0], $
             [0,0,1,1,0,1,0,1,1,0,0], $
             [0,1,1,0,0,1,0,0,1,1,0], $
             [1,0,0,0,0,1,0,0,0,0,1], $
             [0,0,0,0,0,1,0,0,0,0,0], $
             [0,0,0,0,1,1,1,0,0,0,0], $
             [0,0,0,0,0,1,0,0,0,0,0]])
   Im = fltarr(200,200)
   A =  rebin(A,55,50)
   Im[25:79,75:124]=A   
; Rotate through 60 deg, magnify by 1.5 and translate
   Im1 = shift(rot(Im,60,1.5,cubic=-0.5),90,-60)   
; Rotate through 180 deg, shrink by 0.5 and translate
   Im2 = shift(rot(Im,180,0.5,cubic=-0.5),-90,-60)   
; Plot them  
   tvscl,Im+Im1+Im2 
; Invariant moments
   print, hu_moments(Im,/log), format='(7F8.3)'
   print, hu_moments(Im1,/log), format='(7F8.3)'
   print, hu_moments(Im2,/log), format='(7F8.3)'   
   
end