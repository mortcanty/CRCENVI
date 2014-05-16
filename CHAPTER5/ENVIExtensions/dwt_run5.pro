; docformat = 'rst'
; dwt_run5.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

pro dwt_run5_extensions_init
   compile_opt idl2
   e = envi(/current)
   e.addextension, 'DWT Pansharpening', 'dwt_run5', path='Chapter5'
end

;+
; :Description:
;       ENVI extension for panchromatic sharpening
;       under ARSIS model with Mallat's discrete wavelet
;       transform and Daubechies wavelets
; :Params:
;      None                 
; :Uses:
;       ENVI, DWT__DEFINE, ORTHO_REGRESS
; :Author:
;      Mort Canty (2009)         
;-
Pro dwt_run5

COMPILE_OPT idl2

catch, err
   if err ne 0 then begin
      catch, /cancel
      if e5 ne !null then $
         e5.reporterror, 'Error: ' + !error_state.msg $
      else $
         message, !error_state.msg, /continue, /noname
      message, /reset
      return
   endif
   
e5 = envi(/current)
   if e5 eq !null then $
      message, 'This extension requires an interactive ENVI session   
      
print, '---------------------------------'
print, 'DWT Panchromatic sharpening'
print, systime(0)
print, '---------------------------------'
      
; get MS image
   envi_select, title='Select low resolution multi-band input file', $
                       fid=fid, dims=dims, pos=pos
   if (fid eq -1) then begin
     print,'cancelled'
     return
   endif
   print,'Dimensions of low resolution spatial subset'
   print, dims[1:4]
   num_cols = dims[2]-dims[1]+1
   num_rows = dims[4]-dims[3]+1
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
; get number of compressions necessary
   ps1 = map_info1.ps[0]
   num_compr = round(alog(ps/ps1)/alog(2))
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
   if num_compr eq 2 then begin
      X_ul = X_ul-2
      Y_ul = Y_ul-2
   endif

; cutout the corresponding spatial subset of the pan image
   num_cols1 = (2^num_compr)*num_cols
   num_rows1 = (2^num_compr)*num_rows
   dimsp = [-1L,X_ul,X_ul+num_cols1-1,Y_ul,Y_ul+num_rows1-1]
   if (dimsp[1] lt dims1[1]) or (dimsp[3] lt dims1[3]) or $
      (dimsp[2] gt dims1[2]) or (dimsp[4] gt dims1[4]) then $
       message,'PAN and MS do not completely overlap. Try smaller MS'
   pan = envi_get_data(fid=fid1,dims=dimsp,pos=pos1)
   panDWT = Obj_New('DWT',pan)
   panDWT->set_coeff,4
   print,'Dimensions of pan subset'
   print,dimsp[1:4]+1
   str=[file_basename(fname1),strtrim(dimsp[1]+1)+strtrim(dimsp[2]+1)+strtrim(dimsp[3]+1)+strtrim(dimsp[4]+1)]
   envi_info_wid, str, title='Pan spatial subset'

; tie point for fused image
   envi_convert_file_coordinates, fid1, dimsp[1], dimsp[3], ee, n, /to_map
   map_info1.mc[2:3]= [ee,n]

; compress pan image (plus once more for normalization of wavelet coefficients)
   for i=1,num_compr+1 do panDWT->filter
   Xfg = panDWT->get_quadrant(1,/no_edges)
   Xgf = panDWT->get_quadrant(2,/no_edges)
   Xgg = panDWT->get_quadrant(3,/no_edges)
; expand back once
   panDWT->invert
; output arrays for fusion
   case data_type of
    1: image_fused = bytarr(num_cols1,num_rows1,num_bands)
    2: image_fused = intarr(num_cols1,num_rows1,num_bands)
    4: image_fused = fltarr(num_cols1,num_rows1,num_bands)
   12: image_fused = uintarr(num_cols1,num_rows1,num_bands)
   else: begin
         void = dialog_message('format must by byte, (un)signed integer or float',/error)
         return
       end
   endcase
; determine wavelet normalization coefficients
   aa = fltarr(3)
   bb = fltarr(3)
   for k=0,num_bands-1 do begin
;   make a copy of the pan image object and inject the kth MS channel
      msDWT = Obj_New('DWT',panDWT->get_image())
      msDWT->set_coeff,4
      msDWT->set_compressions,num_compr
      ms_band=envi_get_data(fid=fid,dims=dims,pos=pos[k])
      msDWT->inject, ms_band
;   compress once more
      msDWT->filter
;   scatterplots and OLRs of MS WCs vs pan WCs
      Y =  msDWT->get_quadrant(1,/no_edges)
      window,13,xSize=400,ysize=300,xPos=100,yPos= 500,title='FG, Band'+string(k+1)
      plot,Xfg,Y,xRange=[2.0*min(Xfg),2.0*max(Xfg)],yRange=[2.0*min(Y),2.0*max(Y)],pSym=3,background=2^24-1,color=0
      ortho_regress, transpose(Xfg[*]), transpose(Y[*]), a, Xm, Ym, sigma_a, sigma_b
;      aa[0] = a
;      bb[0] = Ym-a*Xm
      oplot,[-100,100],[Ym+a*(-100-Xm),Ym+a*(100-Xm)],pSym=0,color=0
      aa[0] = stddev(Y)/stddev(Xfg)
      bb[0] = mean(Y)-aa[0]*mean(Xfg)
      Y =  msDWT->get_quadrant(2,/no_edges)
      window,14,xSize=400,ysize=300,xPos=100,yPos= 500,title='GF, Band'+string(k+1)
      plot,Xgf,Y,xRange=[2.0*min(Xgf),2.0*max(Xgf)],yRange=[2.0*min(Y),2.0*max(Y)],pSym=3,background=2^24-1,color=0
      ortho_regress, transpose(Xgf[*]), transpose(Y[*]), a, Xm, Ym, sigma_a, sigma_b
;      aa[1] = a
;      bb[1] = Ym-a*Xm
      oplot,[-100,100],[Ym+a*(-100-Xm),Ym+a*(100-Xm)],pSym=0,color=0
      aa[1] = stddev(Y)/stddev(Xgf)
      bb[1] = mean(Y)-aa[1]*mean(Xgf)
      Y =  msDWT->get_quadrant(3,/no_edges)
      window,15,xSize=400,ysize=300,xPos=100,yPos= 500,title='GG, Band'+string(k+1)
      plot,Xgg,Y,xRange=[2.0*min(Xgg),2.0*max(Xgg)],yRange=[2.0*min(Y),2.0*max(Y)],pSym=3,background=2^24-1,color=0
      ortho_regress, transpose(Xgg[*]), transpose(Y[*]), a, Xm, Ym, sigma_a, sigma_b
;      aa[2] = a
;      bb[2] = Ym-a*Xm
      oplot,[-100,100],[Ym+a*(-100-Xm),Ym+a*(100-Xm)],pSym=0,color=0
      aa[2] = stddev(Y)/stddev(Xgg)
      bb[2] = mean(Y)-aa[2]*mean(Xgg)
; expand back
      msDWT->invert
; normalize the MS WCs to the pan WCs
      msDWT->normalize_WC, aa, bb
; expand to full resolution
      for i=1,num_compr do msDWT->invert
; prepare for output
      temp = msDWT->get_image()
      temp = temp > 0.0
      case data_type of
        1: image_fused[*,*,k] = byte(temp < 255.0)
        2: image_fused[*,*,k] = fix(temp < 2047.0)
        4: image_fused[*,*,k] = float(temp)
       12: image_fused[*,*,k] = uint(temp < 2047.0)
      endcase
      Obj_Destroy, msDWT
   endfor ; end loop over bands

   Obj_Destroy, panDWT

; output to memory or file in BSQ
   if (result.outf.in_memory eq 1) then begin
      envi_enter_data, image_fused,map_info=map_info1,xstart=xstart+dims[1], ystart=ystart+dims[3], $
                       bnames='ARSIS_DWT:'+bnames[pos],$
                       descrip='ARSIS_DWT: high-res='+file_basename(fname)+' low-res(ms)='+file_basename(fname1)
      print, 'Result written to memory'
   endif else begin
      openw, unit, result.outf.name, /get_lun
      for i=0,num_bands-1 do writeu, unit, image_fused[*,*,i]
      envi_setup_head ,fname=result.outf.name, ns=num_cols1, $
           nl=num_rows1, nb=num_bands, $
           map_info=map_info1, $
           xstart=xstart+dims[1], ystart=ystart+dims[3], $
           data_type=data_type, interleave=0, /write, $
           bnames='ARSIS_DWT:'+bnames[pos],$
           /open,$
           descrip='ARSIS_DWT: high-res='+file_basename(fname)+' low-res(ms)='+file_basename(fname1)
      print, 'File created ', result.outf.name
      free_lun, unit
   endelse

End

