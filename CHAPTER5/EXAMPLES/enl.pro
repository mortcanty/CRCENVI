pro ENL 
; simple estimation of equivalent
; number of looks (ENL)
   compile_opt idl2
   
; select the polSAR image
   envi_select, title='Enter polSAR image band',$
              /band_only, fid=fid, dims=dims, pos=pos
   if (fid eq -1) then begin
      print,'cancelled'
      return
   endif
   envi_file_query, fid, fname=fname
   Gs = envi_get_data(/complex, dims=dims, $
                      fid=fid, pos=pos)  
                        
; get the ROI pixels 
   roi_ids = envi_get_roi_ids(fid=fid)
   idx = envi_get_roi(roi_ids[0])
   
; calculate ENL                             
   mn = mean(abs(Gs[idx]))
   var = variance(abs(Gs[idx]))
   print, 'Input file: ' + fname
   print, 'ENL = ', real_part(mn^2/var) 
 
end