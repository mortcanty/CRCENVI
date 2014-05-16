; docformat = 'rst'
; enltm_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
    
PRO enltm_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'ENL TM estimation', $
      REF_VALUE = 'AIRSAR Scattering Classification', $
      EVENT_PRO = 'enltm_run', $
      UVALUE = 'ENL',$
      POSITION = 'after', $
      /SEPARATOR
END

pro get_windex, j, cols, windex
    windex = lonarr(49)
    windex[0:6]   = (j-3)*cols + [0,1,2,3,4,5,6]
    windex[7:13]  = (j-2)*cols + [0,1,2,3,4,5,6]
    windex[14:20] = (j-1)*cols + [0,1,2,3,4,5,6]
    windex[21:27] = (j)*cols   + [0,1,2,3,4,5,6] 
    windex[28:34] = (j+1)*cols + [0,1,2,3,4,5,6] 
    windex[35:41] = (j+2)*cols + [0,1,2,3,4,5,6]
    windex[42:48] = (j+3)*cols + [0,1,2,3,4,5,6]   
end


;+
; :Description:   
;    Estimation of ENL for polSAR covariance images 
;    with trace moment method
;    Uses diagonal elements only
;    Anfinsen et al. (2009) IEEE TGARS 47(11), 3795-3809
; :Params:
;      event:  in, required 
;         if called from the ENVI Classic menu   
; :KEYWORDS:
;    NONE    
; :Uses:
;    POLSAR, ENVI    
; :Author:
;       Mort Canty (2013)        
;-
pro enltm_run, event 

   Compile_Opt idl2
;  Standard error handling.
   Catch, theError
   IF theError NE 0 THEN BEGIN
      Catch, /CANCEL
      void = Error_Message()
      RETURN
   ENDIF    
    
   print, '--------------------------------------------'
   print, 'ENL TM estimation for polSAR covariance image'
   print, systime(0)
   print, '--------------------------------------------'        
; input multi-look averaged covariance image  
   fname1 = envi_pickfile(filter='*.hdr',title='Select SAR image header')
   if fname1 eq '' then return
   base = widget_auto_base(title='Byte order')
   list = ['Little Endian', 'Big Endian']
   wm = widget_menu(base,list=list,default_ptr=0,uvalue='menu',/excl,/auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then return
   if result.menu eq 0 then endian = 'little' else endian = 'big'    
   ps = POLSAR(fname1,endian=endian)         
   ps.getProperty, lines=rows, samples=cols, $
     c11=c11, c22=c22, c33=c33, map_info=map_info
   C = fltarr(3,cols*rows)
   C[0,*] = real_part(c11)
   C[1,*] = real_part(c22)
   C[2,*] = real_part(c33)  
   ENL_TM = fltarr(cols,rows)
   start = systime(2) 
   iter = 0
   print, 'polSAR image:  '+fname1
   progressbar = Cgprogressbar(/cancel,title='ENL...')  
   progressbar->start   
   for j = 3L,rows-4 do begin
      if progressbar->CheckCancel() then begin
         progressbar->destroy
         obj_destroy, ps
         return
      endif 
      get_windex, j, cols, windex      
      for i = 3L,cols-4 do begin  
         wind = C[*,windex]
;       check for zero data
         if total(wind) gt 0 then begin
;          calculate the TM (trace moment) estimate          
            EC = total(wind,2)/49       
            TrECEC = total(EC*EC)           
            ETrCsqr = total( total(wind,1)^2 )/49           
            TrECsqr = total(EC)^2
            ENL_TM[i,j] = TrECEC/(ETrCsqr - TrECsqr)      
         endif      
         windex++  
      endfor 
      progressbar->Update,j*100/rows  
   endfor   
   progressbar->destroy 
   idx = where(ENL_TM gt 0)  
   print, 'Iteration:  ',strtrim(iter+1,2), $
     '   Mean TM: ', mean(ENL_TM[idx])
   obj_destroy, ps 
   p1 = plot(histogram(ENL_TM[*],min=0,max=50,nbins=50))
   envi_enter_data, ENL_TM, bnames = ['ENL_TM'], map_info = map_info  
   print, 'done, elapsed time: '+strtrim(systime(2)-start,2)+' sec'    
end