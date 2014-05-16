; docformat = 'rst'
; atwt_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO atwt_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'A Trous Wavelet Transform', $
      REF_VALUE = 'CN Spectral Sharpening', $
      EVENT_PRO = 'atwt_run', $
      UVALUE = 'ATWT',$
      POSITION = 'after'
END

;+
; :Description:
;       ENVI extension for panchromatic sharpening under
;       ARSIS model with "a trous" wavelet transform.
; :Params:
;      event:  in, required 
;         if called from the ENVI menu              
; :Uses:
;       ENVI, ATWT__DEFINE, ORTHO_REGRESS
; :Author:
;      Mort Canty (2013)             
;-
Pro atwt_run, event

COMPILE_OPT IDL2

   print, '-----------------------------------'
   print, 'ARSIS ATWT Model for Wavelet Fusion'
   print, systime(0)
   print, '-----------------------------------'

catch, theError
if theError ne 0 then begin
   catch,/cancel
   ok = dialog_message(!Error_State.msg,/error)
   return
endif

; get MS image
   envi_select, title='Select  low resolution multi-band input file', $
                       fid=fid, dims=dims, pos=pos
   if (fid eq -1) then begin
     print,'cancelled'
     return
   endif
   print,'Dimensions of low resolution spatial subset'
   print, dims[1:4]
   num_cols = dims[2]-dims[1]+1L
   num_rows = dims[4]-dims[3]+1L
   num_bands = n_elements(pos)
   envi_file_query, fid, fname=fname,bnames=bnames,data_type=data_type,xstart=xstart,ystart=ystart
   map_info = envi_get_map_info(fid=fid)
   ps = map_info.ps[0]
   print,'Low resolution (multispectral) image chosen: ',fname
   print,'Bands ',pos+1

; get pan image
   envi_select, title='Select high resolution (pan) input band',fid=fid1,pos=pos1,dims=dims1,/band
   if (fid1 eq -1) then begin
     print,'cancelled'
     return
   endif
   envi_file_query, fid1, fname=fname1,data_type=data_type1,xstart=xstart1,ystart=ystart1
   map_info1 = envi_get_map_info(fid=fid1)
; get number of applications necessary
   ps1 = map_info1.ps[0]
   num_trans=round(alog(ps/ps1)/alog(2))
   ratio = round(ps/ps1)
   print,'High resolution image chosen: ',fname1
   print,'Band ',pos1+1

; get output destination
   base = widget_auto_base(title='ARSIS Fusion Output')
   sb = widget_base(base, /row, /frame)
   wp = widget_outfm(sb, uvalue='outf', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
     print, 'cancelled'
     return
   endif

   widget_control,/hourglass

; find upper left position of low-res image within pan image
   envi_convert_file_coordinates,fid,dims[1],dims[3],e,n,/to_map
   envi_convert_file_coordinates,fid1,X_ul,Y_ul,e,n
   X_ul = round(X_ul)
   Y_ul = round(Y_ul)
; shift to pixel center
   if num_trans eq 2 then begin
      X_ul = X_ul;-2
      Y_ul = Y_ul;-2
   endif

; cutout the corresponding spatial subset of the pan image
   num_cols1 = (2^num_trans)*num_cols
   num_rows1 = (2^num_trans)*num_rows

   dimsp = [-1L,X_ul,X_ul+num_cols1-1,Y_ul,Y_ul+num_rows1-1]
   if (dimsp[1] lt dims1[1]) or (dimsp[3] lt dims1[3]) or $
      (dimsp[2] gt dims1[2]) or (dimsp[4] gt dims1[4]) then $
       message,'PAN and MS do not completely overlap. Try smaller MS'      
   pan = envi_get_data(fid=fid1,dims=dimsp,pos=pos1)
   print,'Dimensions of pan subset'
   print,dimsp[1:4]+1
   str=[file_basename(fname1),strtrim(dimsp[1]+1)+strtrim(dimsp[2]+1)+strtrim(dimsp[3]+1)+strtrim(dimsp[4]+1)]
   envi_info_wid, str, title='Pan spatial subset'

; tie point for fused image
   envi_convert_file_coordinates, fid1, dimsp[1], dimsp[3], e, n, /to_map
   map_info1.mc[2:3]= [e,n]

; output arrays for fusion
   case data_type of
    1: image_fused = bytarr(num_cols1,num_rows1,num_bands)
    2: image_fused = intarr(num_cols1,num_rows1,num_bands)
    4: image_fused = fltarr(num_cols1,num_rows1,num_bands)
   12: image_fused = uintarr(num_cols1,num_rows1,num_bands)
   else: begin
         void = dialog_message('format must be byte, float or (un)signed integer',/error)
         return
       end
   endcase
   smpl = randomu(seed,100000L,/long) mod (num_cols1*num_rows1)

; loop over MS bands
   for k=0,num_bands-1 do begin
; transform the pan image
      msATWT = Obj_New('ATWT',pan)
      for i=1,num_trans do msATWT->compress
; sample details for pan and then transform back to save space
      X=(msATWT->get_image(num_trans))[smpl]
      for i=1,num_trans do msATWT->expand
; resize the ms band to scale of the pan image using NN resampling
      ms_band=envi_get_data(fid=fid,dims=dims,pos=pos[k])
      envi_enter_data, ms_band,r_fid=r_fid
      dims2 = [-1L,0,num_cols-1,0,num_rows-1]
      envi_doit, 'resize_doit',/in_memory,fid=r_fid,pos=0,dims=dims2,r_fid=r_fid1,interp=0,rfact=[1.0,1.0]/2^num_trans
      dims3 = [-1L,0,num_cols1-1,0,num_rows1-1]
      ms_band=envi_get_data(fid=r_fid1,dims=dims3,pos=0)
      envi_file_mng, id=r_fid,/remove
      envi_file_mng, id=r_fid1,/remove
; sample details for ms
      tempATWT = Obj_New('ATWT',ms_band)
      for i=1,num_trans do tempATWT->compress
      Y=(tempATWT->get_image(num_trans))[smpl]
; get image for injection
      img = tempATWT->get_image(0)
      Obj_destroy, tempATWT
; scatterplots and OLRs of ms details vs pan details
      window,13,xSize=400,ysize=300,xPos=100,yPos= 500,title='Band '+strtrim(k+1,2)
      plot,X,Y,xRange=[2.0*min(X),2.0*max(X)],yRange=[2.0*min(Y),2.0*max(Y)],pSym=3,background=2^24-1,color=0
;      ortho_regress, transpose(X), transpose(Y), a, Xm, Ym, sigma_a, sigma_b
;      oplot,[-200,200],[Ym+a*(-200-Xm),Ym+a*(200-Xm)],pSym=0,color=0
      aa = stddev(Y)/stddev(X)
      bb = mean(Y)-aa*mean(X)
; transform again and inject the filtered ms band
      for i=1,num_trans do msATWT->compress
      msATWT->inject, img
; normalize the details
      msATWT->normalize_wc,aa,bb
; expand to full resolution
      for i=1,num_trans do msATWT->expand
; prepare for output
      temp = msATWT->get_image(0)
      temp = temp > 0.0
      case data_type of
        1: image_fused[*,*,k] = byte(temp < 255.0)
        2: image_fused[*,*,k] = fix(temp < 2047.0)
        4: image_fused[*,*,k] = float(temp)
       12: image_fused[*,*,k] = uint(temp < 2047.0)
      endcase
      Obj_destroy, msATWT
   endfor ; end loop over bands

; output to memory or file in BSQ
   if (result.outf.in_memory eq 1) then begin
      envi_enter_data, image_fused,map_info=map_info1,xstart=xstart+dims[1], ystart=ystart+dims[3], $
                       bnames='ARSIS_ATWT:'+bnames[pos],$
                       descrip='ARSIS_ATWT: high-res='+file_basename(fname)+' low-res(ms)='+file_basename(fname1)
      print, 'Result written to memory'
   endif else begin
      openw, unit, result.outf.name, /get_lun
      for i=0,num_bands-1 do writeu, unit, image_fused[*,*,i]
      envi_setup_head ,fname=result.outf.name, ns=num_cols1, $
           nl=num_rows1, nb=num_bands, $
           map_info=map_info1, $
           xstart=xstart+dims[1], ystart=ystart+dims[3], $
           data_type=data_type, interleave=0, /write, $
           bnames='ARSIS_ATWT:'+bnames[pos],$
           /open,$
           descrip='ARSIS_ATWT: high-res='+file_basename(fname)+' low-res(ms)='+file_basename(fname1)
      print, 'File created ', result.outf.name
      free_lun, unit
   endelse

End

