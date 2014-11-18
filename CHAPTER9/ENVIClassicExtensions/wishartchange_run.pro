; docformat = 'rst'
; wishartchange_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
    
PRO wishartchange_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Complex Wishart Change Detection', $
      REF_VALUE = 'AIRSAR Scattering Classification', $
      EVENT_PRO = 'wishartchange_run', $
      UVALUE = 'WISHART',$
      POSITION = 'after'
END
    
;+
; :Description:
;    Polarimetric change detection for polarised SAR images
;    (multi-look covariance representation)  
;    Perform test for equality of two complex 1x1, 2x2 or 3x3
;    variance-covariance matrices representing
;      full polarimetry (9 bands)
;      dual polarimetry (4 bands
;      single polarimetry (1 band)
;    Based on a Matlab script by Allan Nielsen.      
; :Params:
;    event: in, optional
;       required if called from ENVI           
; :KEYWORDS:
;    None    
; :Uses:
;    ENVI  
; :Author:
;    Mort Canty (2014)        
;-
pro wishartchange_run, event 

   Compile_Opt idl2
;  Standard error handling.
   Catch, theError
   IF theError NE 0 THEN BEGIN
      Catch, /CANCEL
      void = Error_Message()
      RETURN
   ENDIF    
    
   print, '-------------------------------------------------------------'
   print, 'Complex Wishart change detection for polSAR covariance images'
   print, systime(0)
   print, '-------------------------------------------------------------'        
; input multi-look averaged covariance images and ENLs
   envi_select, title='Choose first covariance matrix image', $
     fid=fid1, dims=dims1, pos=pos1, /no_spec, /no_dims
   if (fid1 eq -1) then begin
     print, 'cancelled'
     return
   end
   envi_file_query, fid1, fname=fname1, bnames=bnames1
   print, 'First image: '+fname1 
   cols1 = dims1[2]-dims1[1]+1
   rows1 = dims1[4]-dims1[3]+1
   bands1 = n_elements(pos1)
; bands1 = 9: quad pol
; bands1 = 4  dual pol
; bands1 = 1  single pol
   base = widget_auto_base(title='Number of looks, first image')
   we = widget_param(base, dt=4, field=5, floor=1, $
     default=8, uvalue='param', xsize=32, /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then n1=float(8) else n1=float(result.param)
   print, 'Number of looks: '+strtrim(n1,2)
   envi_select, title='Choose second covariance matrix image', $
     fid=fid2, dims=dims2, pos=pos2, /no_spec, /no_dims
   if (fid2 eq -1) then begin
     print, 'cancelled'
     return
   end
   envi_file_query, fid2, fname=fname2, bnames=bnames2
   print, 'First image: '+fname2 
   cols2 = dims2[2]-dims2[1]+1
   rows2 = dims2[4]-dims2[3]+1
   bands2 = n_elements(pos2)
   base = widget_auto_base(title='Number of looks, second image')
   we = widget_param(base, dt=4, field=5, floor=1, $
     default=n1, uvalue='param', xsize=32, /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then n2=float(8) else n2=float(result.param) 
   print, 'Number of looks: '+strtrim(n2,2)
   base = widget_auto_base(title='Significance level')
   we = widget_param(base, dt=4, field=5, floor=0., $
   default=0.01, uvalue='param', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then siglevel=0.01 else siglevel = float(result.param) 
   print, 'Change significance level:' + strtrim(siglevel,2)  
     
   base = widget_auto_base(title='Output image')
   sb = widget_base(base, /row, /frame)
   wp = widget_outfm(sb, uvalue='outf', /auto)
   result = auto_wid_mng(base)
   
   if not result.accept then begin
      print, 'Output was cancelled'
      return
   endif   
   
   map_info = envi_get_map_info(fid=fid1) 
   envi_convert_file_coordinates, fid1, $
     dims1[1], dims1[3], e, n, /to_map
   map_info.mc = [0D,0D,e,n]
   
; get tie points for geographic matching
   pts = fltarr(4,4)
   x1 = dims1[1]
   y1 = dims1[3]
   envi_convert_file_coordinates, fid1, x1, y1, e, n, /to_map
   envi_convert_file_coordinates, fid2, x2, y2, e, n
   pts[*,0] = [x1,y1,x2,y2]
   x1 = dims1[2]
   y1 = dims1[3]
   envi_convert_file_coordinates, fid1, x1, y1, e, n, /to_map
   envi_convert_file_coordinates, fid2, x2, y2, e, n
   pts[*,1] = [x1,y1,x2,y2]
   x1 = dims1[1]
   y1 = dims1[4]
   envi_convert_file_coordinates, fid1, x1, y1, e, n, /to_map
   envi_convert_file_coordinates, fid2, x2, y2, e, n
   pts[*,2] = [x1,y1,x2,y2]
   x1 = dims1[2]
   y1 = dims1[4]
   envi_convert_file_coordinates, fid1, x1, y1, e, n, /to_map
   envi_convert_file_coordinates, fid2, x2, y2, e, n
   pts[*,3] = [x1,y1,x2,y2]
   
 ; co-register  
 ;-----workaround------------------
   x0 = map_info.mc[2] 
   y0 = map_info.mc[3]
   xsize = map_info.ps[0]*cols1 
   ysize = map_info.ps[1]*rows1 
 ;---------------------------------   
   ENVI_DOIT, 'ENVI_REGISTER_DOIT', b_fid=fid1, w_fid=fid2, $
                                    w_dims=dims2, $
                                    w_pos=pos2, $
                                    pts=pts, r_fid=r_fid, $
                                    method=0, $
                                    x0 = x0, y0 = y0, $
                                    xsize = xsize, ysize = ysize, $
                                    out_name = fname2+'_reg'                                   
   fid2 = r_fid  
   dims2 = [-1,0,cols1-1,0,rows1-1]                                   
; calculate determinants and change statistic distribution parameters  
   case bands1 of
     9: begin
       k1 = n1*envi_get_data(fid=fid1,dims=dims1,pos=0)  ;c11
       a1 = envi_get_data(fid=fid1,dims=dims1,pos=1)  ;c12
       im = envi_get_data(fid=fid1,dims=dims1,pos=2)
       a1 = n1*complex(a1,im)
       rho1 = envi_get_data(fid=fid1,dims=dims1,pos=3) ;c13
       im = envi_get_data(fid=fid1,dims=dims1,pos=4)
       rho1 = n1*complex(rho1,im)
       xsi1 = n1*envi_get_data(fid=fid1,dims=dims1,pos=5) ;c22
       b1 = envi_get_data(fid=fid1,dims=dims1,pos=6)   ;c23
       im = envi_get_data(fid=fid1,dims=dims1,pos=7)
       b1 = n1*complex(b1,im)
       zeta1 = n1*envi_get_data(fid=fid1,dims=dims1,pos=8);c33
       det1 = k1*xsi1*zeta1 + 2*real_part(a1*b1*conj(rho1)) - xsi1*(abs(rho1)^2) - k1*(abs(b1)^2) - zeta1*(abs(a1)^2)
       span1 = k1 + 2*xsi1 + zeta1
       k2 = n2*envi_get_data(fid=fid2,dims=dims2,pos=0)  ;c11
       a2 = envi_get_data(fid=fid2,dims=dims2,pos=1)  ;c12
       im = envi_get_data(fid=fid2,dims=dims2,pos=2)
       a2 = n2*complex(a2,im)
       rho2 = envi_get_data(fid=fid2,dims=dims2,pos=3) ;c13
       im = envi_get_data(fid=fid2,dims=dims2,pos=4)
       rho2 = n2*complex(rho2,im)
       xsi2 = n2*envi_get_data(fid=fid2,dims=dims2,pos=5) ;c22
       b2 = envi_get_data(fid=fid2,dims=dims2,pos=6)   ;c23
       im = envi_get_data(fid=fid2,dims=dims2,pos=7)
       b2 = n2*complex(b2,im)
       zeta2 = n2*envi_get_data(fid=fid2,dims=dims2,pos=8);c33
       det2 = k2*xsi2*zeta2 + 2*real_part(a2*b2*conj(rho2)) - xsi2*(abs(rho2)^2) - k2*(abs(b2)^2) - zeta2*(abs(a2)^2)       
       k3    = k1 + k2
       a3    = a1 + a2
       rho3  = rho1 +  rho2
       xsi3  = xsi1 +  xsi2
       b3    = b1 +    b2
       zeta3 = zeta1 + zeta2  
       det3 = k3*xsi3*zeta3 + 2*real_part(a3*b3*conj(rho3)) - xsi3*(abs(rho3)^2) - k3*(abs(b3)^2) - zeta3*(abs(a3)^2)  
       p = 3
       f = p^2
       cst = p*((n1+n2)*alog(n1+n2)-n1*alog(n1)-n2*alog(n2))
       rho = 1. - (2.*p^2-1.)*(1./n1 + 1./n2 - 1./(n1+n2))/(6.*p)
       omega2 = -(p*p/4.)*(1. - 1./rho)^2 + p^2*(p^2-1.)*(1./n1^2 + 1./n2^2 - 1./(n1+n2)^2)/(24.*rho^2)   
     end
     4: begin 
       k1 = n1*envi_get_data(fid=fid1,dims=dims1,pos=0) ;c11
       a1 = envi_get_data(fid=fid1,dims=dims1,pos=1) ;c12
       im = envi_get_data(fid=fid1,dims=dims1,pos=2)
       a1 = n1*complex(a1,im)
       xsi1 = n1*envi_get_data(fid=fid1,dims=dims1,pos=3);c22
       det1 = k1*xsi1 - abs(a1)^2
       span1 = k1 + xsi1
       k2 = n2*envi_get_data(fid=fid2,dims=dims2,pos=0) ;c11
       a2 = envi_get_data(fid=fid2,dims=dims2,pos=1) ;c12
       im = envi_get_data(fid=fid2,dims=dims2,pos=2)
       a2 = n2*complex(a2,im)
       xsi2 = n2*envi_get_data(fid=fid2,dims=dims2,pos=3);c22
       det2 = k2*xsi2 - abs(a2)^2       
       k3    = k1 + k2
       a3    = a1 + a2
       xsi3  = xsi1 +  xsi2
       det3 = k3*xsi3 - abs(a3)^2  
       p = 2
       cst = p*((n1+n2)*alog(n1+n2)-n1*alog(n1)-n2*alog(n2))
       f = p^2;
       rho = 1-(2*f-1)*(1./n1+1./n2-1./(n1+n2))/(6.*p);
       omega2 = -f/4.*(1-1./rho)^2 + f*(f-1)*(1./n1^2+1./n2^2-1./(n1+n2)^2)/(24.*rho^2);
     end
     1: begin
       k1 = n1*envi_get_data(fid=fid1,dims=dims1,pos=0) ;c11
       span1 = k1
       det1 = k1
       k2 = n2*envi_get_data(fid=fid2,dims=dims2,pos=0) ;c11
       det2 = k2
       k3 = k1 + k2
       det3 = k3
       p = 1
       cst = p*((n1+n2)*alog(n1+n2)-n1*alog(n1)-n2*alog(n2))
       f = p^2;
       rho = 1-(2.*f-1)*(1./n1+1./n2-1./(n1+n2))/(6.*p);
       omega2 = -f/4.*(1-1./rho)^2+f*(f-1)*(1./n1^2+1./n2^2-1./(n1+n2)^2)/(24.*rho^2);
     end
     else: begin
      message, 'Incorrect input data, aborting...'
      return  
     end          
   endcase
   envi_file_mng, id = fid2, /remove
; chck fo bad data   
   idx = where(det1 le 0,count)
   if count gt 0 then begin
      print, 'Warning: det(image1) has non-positive values'
      det1[idx] = 0.0001;(machar()).eps
   endif   
   idx = where(det2 le 0,count)
   if count gt 0 then begin
      print, 'Warning: det(image2) has non-positive values'
      det2[idx] = 0.0001;(machar()).eps
   endif           
   idx = where(det3 le 0,count)
   if count gt 0 then begin
      print, 'Warning: det(image1+image2) has non-positive values'
      det3[idx] = 0.0001;(machar()).eps
   endif   
; change statistic             
   lnQ = cst+n1*alog(det1)+n2*alog(det2)-(n1+n2)*alog(det3)
   Z = -2*rho*lnQ
; change probability   
   PP =  (1.-omega2)*chisqr_pdf(Z,f)+omega2*chisqr_pdf(Z,f+4)  
; output 
   outimage = fltarr(cols1,rows1,4)
   outimage[*,*,0] = span1
   outimage[*,*,1] = Z
   outimage[*,*,2] = PP
   PP = median(PP,3)
   idx = where(PP lt 1.0-siglevel)
   PP[idx] = 0.0
   outimage[*,*,3] = PP
; median filter on change map-----   
   outimage[*,*,3] = median(PP,3)
; --------------------------------   
   bnames = ['span t1', '-2*rho*lnQ', 'change prob', 'change ('+strtrim(siglevel,2)+' significance level)'] 
; change map as ENVI classification image   
   classimage = byte(outimage[*,*,3]+0.1)  
   envi_enter_data, classimage,$
      file_type = 3, $
      num_classes=2, $
      map_info=map_info, $
      bnames = ['change ('+strtrim(siglevel,2)+' significance level)'], $
      class_names=['no change','change ('+strtrim(siglevel,2)+' significancee level)'], $
      lookup = [[0,0,0],[255,0,0]]
; write to memory or disk  
   if (result.outf.in_memory eq 1) then begin
      envi_enter_data, outimage, bnames=bnames, map_info=map_info, $
         descrip='Wishart change( '+file_basename(fname1)+','+file_basename(fname2)+' )'
      print, 'Change image written to memory'
   end else begin
      openw, unit, result.outf.name, /get_lun
      writeu, unit, outimage
      envi_setup_head,fname=result.outf.name, ns=cols1, nl=rows1, nb=4, $
                    data_type=4, $
                    interleave=0, $
                    file_type=0, $  
                    bnames = bnames, $  
                    map_inf=map_info, $
                    descrip='Wishart change( '+file_basename(fname1)+','+file_basename(fname2)+' )', $
                    /write,/open
      print, 'File created ', result.outf.name
      free_lun, unit
   endelse
   
end