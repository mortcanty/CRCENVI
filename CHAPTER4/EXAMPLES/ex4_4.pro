pro EX4_4

envi_select, title='Choose multispectral band', $
             fid=fid, dims=dims,pos=pos, /band_only
if (fid eq -1) then begin
   print, 'cancelled'
   return
endif

; get the image band from ENVI
g = envi_get_data(fid=fid,dims=dims,pos=pos[0])
; create an instance of the DWT class with image band
aDWT = Obj_New('DWT', g)
; apply the filter bank once
aDWT->filter
; return the result to ENVI
envi_enter_data, aDWT->get_image()
; apply again
aDWT->filter
; return the result to ENVI
envi_enter_data, aDWT->get_image()
; remove the class instance from the heap
Obj_Destroy, aDWT

end

