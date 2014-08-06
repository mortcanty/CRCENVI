; docformat = 'rst'
; enlml_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
    
PRO enlml_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'ENL estimation', $
      REF_VALUE = 'AIRSAR Scattering Classification', $
      EVENT_PRO = 'enl_ml_run', $
      UVALUE = 'ENLML',$
      POSITION = 'after'
END

pro get_windex3, j, cols, windex
    windex = lonarr(9)
    windex[0:2] = (j-1)*cols + [0,1,2]
    windex[3:5] = (j  )*cols + [0,1,2]
    windex[6:8] = (j+1)*cols + [0,1,2]
end

pro get_windex7, j, cols, windex
    windex = lonarr(49)
    windex[0:6]   = (j-3)*cols + [0,1,2,3,4,5,6]
    windex[7:13]  = (j-2)*cols + [0,1,2,3,4,5,6]
    windex[14:20] = (j-1)*cols + [0,1,2,3,4,5,6]
    windex[21:27] = (j  )*cols + [0,1,2,3,4,5,6] 
    windex[28:34] = (j+1)*cols + [0,1,2,3,4,5,6] 
    windex[35:41] = (j+2)*cols + [0,1,2,3,4,5,6]
    windex[42:48] = (j+3)*cols + [0,1,2,3,4,5,6]   
end

pro get_windex11, j, cols, windex
    windex = lonarr(121)
    windex[0:10]  =   (j-5)*cols + [0,1,2,3,4,5,6,7,8,9,10]
    windex[11:21] =   (j-4)*cols + [0,1,2,3,4,5,6,7,8,9,10]
    windex[22:32] =   (j-3)*cols + [0,1,2,3,4,5,6,7,8,9,10]
    windex[33:43] =   (j-2)*cols + [0,1,2,3,4,5,6,7,8,9,10]
    windex[44:54] =   (j-1)*cols + [0,1,2,3,4,5,6,7,8,9,10]    
    windex[55:65] =   (j  )*cols + [0,1,2,3,4,5,6,7,8,9,10]
    windex[66:76] =   (j+1)*cols + [0,1,2,3,4,5,6,7,8,9,10]
    windex[77:87] =   (j+2)*cols + [0,1,2,3,4,5,6,7,8,9,10]
    windex[88:98] =   (j+3)*cols + [0,1,2,3,4,5,6,7,8,9,10]
    windex[99:109] =  (j+4)*cols + [0,1,2,3,4,5,6,7,8,9,10]
    windex[110:120] = (j+5)*cols + [0,1,2,3,4,5,6,7,8,9,10]         
end    
    
;+
; Return digamma function (derivative of gamma function divided by
; gamma function)
; :Params:
;    x, in, required
;       Argument for digamma function
;    eps, in, optional, default=1e-12
;       Precision for result
; :History:
;    15 Apr 2008 Written, Anthony Smith
;-
FUNCTION ajs_digamma, x, eps
  compile_opt idl2

  x = double(x)
  IF n_elements(eps) EQ 0 THEN $
     eps = 1e-12
  
  ;; Euler-Mascheroni constant
  gamma = 0.57721566490153286060651209008240243104215933593992d

  psi = - gamma
  n = 1
  REPEAT BEGIN
      delta_psi = (x - 1) / (n * (n + x - 1))
      psi += delta_psi
      n += 1
  ENDREP UNTIL (abs(delta_psi) LT eps)

  return, psi
END

;+
; Construct lookup table for ML ENL estimation
;-
pro Make_Lookup
  compile_opt idl2
   Llu = fltarr(3,800)
   digamma = fltarr(800)
   for L = 10,799 do begin
       digamma[L] = ajs_digamma(L/10.0+0.0)
       if (L mod 10) eq 0 then print, L/10
   end    
   for L = 10,799 do Llu[0,L] =  -digamma[L] + alog(L/10.0)
   for L = 20,799 do Llu[1,L] =  -(2*digamma[L-10] + 1.0/(L/10.0-1.0)) + 2*alog(L/10.0)
   for L = 30,799 do Llu[2,L] =  -(3*digamma[L-20] + 2.0/(L/10.0-2.0) + 1.0/(L/10.0-1.0)) + 3*alog(L/10.0)
   p = plot(LLu[0,*])
   p = plot(LLu[1,*],/overplot)
   p = plot(LLu[2,*],/overplot)
   openw, lun, 'd:\idl\crc\lookup.txt', /get_lun
   printf, lun, LLu
   free_lun, lun
end
    
;+
; :Description:   
;    Estimation of ENL for polSAR covariance images
;    using ML method with full covariance matrix (quad, dual or single)
;    Anfinsen et al. (2009) IEEE TGARS 47(11), 3795-3809
;    Takes input from covariance matrix format images generated
;    from polsaringest_run.pro
; :Params:
;      event:  in, required 
;         if called from the ENVI Classic menu   
; :KEYWORDS:
;    NONE    
; :Uses:
;    COYOTE, ENVI    
; :Author:
;       Mort Canty (2014)        
;-
pro enlml_run, event 

   Compile_Opt idl2
;  Standard error handling.
   Catch, theError
   IF theError NE 0 THEN BEGIN
      Catch, /CANCEL
      void = Error_Message()
      RETURN
   ENDIF    
    
   print, '------------------------------------------'
   print, 'ENL estimation for polSAR covariance image'
   print, systime(0)
   print, '------------------------------------------'        
; input multi-look averaged covariance image  
   envi_select, title='Choose (spatial subset of) covariance matrix image', $
     fid=fid, dims=dims, pos=pos
   if (fid eq -1) then begin
     print, 'cancelled'
     return
   end
   envi_file_query, fid, fname=fname, bnames=bnames
   cols = dims[2]-dims[1]+1
   rows = dims[4]-dims[3]+1
   bands = n_elements(pos)
; bands = 9: quad pol
; bands = 4  dual pol
; bands = 1  single pol
   
   base = widget_auto_base(title='Window size')
   list = ['3x3', '7x7', '11x11']
   wm = widget_menu(base, list=list, uvalue='menu', default_ptr=1, /excl, /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then return
   case result.menu of
      0: ws = 3
      1: ws = 7
      2: ws = 11
   endcase  
   
   print, 'polSAR image:  '+fname
   
   case bands of 
      9: begin
           print, bnames
           k = envi_get_data(fid=fid,dims=dims,pos=0)  ;c11
           a = envi_get_data(fid=fid,dims=dims,pos=1)  ;c12
           im = envi_get_data(fid=fid,dims=dims,pos=2) 
           a = complex(a,im)
           rho = envi_get_data(fid=fid,dims=dims,pos=3) ;c13
           im = envi_get_data(fid=fid,dims=dims,pos=4)
           rho = complex(rho,im)
           xsi = envi_get_data(fid=fid,dims=dims,pos=5) ;c22          
           b = envi_get_data(fid=fid,dims=dims,pos=6)   ;c23
           im = envi_get_data(fid=fid,dims=dims,pos=7)
           b = complex(b,im)   
           zeta = envi_get_data(fid=fid,dims=dims,pos=8);c33  
           det = k*xsi*zeta + 2*real_part(a*b*conj(rho)) - xsi*(abs(rho)^2) - k*(abs(b)^2) - zeta*(abs(a)^2)
           d = 3
        end
      4: begin
           print, 'c11-> '+bnames[0]+' c12re-> '+bnames[1]+' c12im-> '+bnames[2]+' c22-> '+bnames[3]
           k = envi_get_data(fid=fid,dims=dims,pos=0) ;c11
           a = envi_get_data(fid=fid,dims=dims,pos=1) ;c12
           im = envi_get_data(fid=fid,dims=dims,pos=2)
           a = complex(a,im)
           xsi = envi_get_data(fid=fid,dims=dims,pos=3);c22
           det = k*xsi - abs(a)^2
           d = 2
         end  
      1: begin
           print, 'c11-> '+bnames[0]
           k = envi_get_data(fid=fid,dims=dims,pos=0) ;c11 
           det = k
           d = 1
         end   
   endcase
    
   case ws of
      3: st = 1L
      7: st = 3L
      11: st = 5L
   endcase   
 
   ENL_ML = fltarr(cols,rows)
   start = systime(2) 
   lu = file_which('lookup.txt',/include_current_dir)
   openr, lun, lu, /get_lun
   LLu = fltarr(3,800)
   readf, lun, LLu
   free_lun, lun
   progressbar = Obj_New('cgprogressbar', /cancel, title='ENL ...') 
   progressbar->start   
   for j = st,rows-st-1 do begin
      if progressbar->CheckCancel() then begin
         progressbar->destroy
         return
      endif 
      case ws of 
         3: get_windex3, j, cols, windex    
         7: get_windex7, j, cols, windex
         11: get_windex11, j, cols, windex
         endcase       
      for i = st,cols-st-1 do begin 
         detC = det[windex]
;       check for zero data
         if min(detC) gt 0 then begin
;          calculate the ML (maximum likelihood) estimate
            avlogdetC = total(alog(detC))/ws^2 
            avdetC = total(detC)/ws^2 
            case bands of
               9: begin
                    k1 = total(k[windex])/ws^2
                    a1 = total(a[windex])/ws^2
                    rho1 = total(rho[windex])/ws^2
                    xsi1 = total(xsi[windex])/ws^2
                    b1 = total(b[windex])/ws^2
                    zeta1 = total(zeta[windex])/ws^2
                    detavC = k1*xsi1*zeta1 + 2*real_part(a1*b1*conj(rho1)) - xsi1*(abs(rho1)^2) - k1*(abs(b1)^2) - zeta1*(abs(a1)^2)
                  end 
               4: begin
                    k1 = total(k[windex])/ws^2
                    xsi1 = total(xsi[windex])/ws^2
                    a1 = total(a[windex])/ws^2
                    detavC = k1*xsi1 - abs(a1)^2          
                  end  
               1: detavC = total(k[windex])/ws^2                                              
            endcase
            logdetavC = alog(detavC)
            arr =  avlogdetC - logdetavC + Llu[d-1,*]
            ell = where(arr*shift(arr,1)lt 0,count) 
            if count gt 1 then ENL_ML[i,j] = float(ell[-1])/10.0 $
               else ENL_ML[i,j] = ENL_ML[(i-1)>0,j]            
         endif      
         windex++  
      endfor 
      progressbar->Update,j*100/rows  
   endfor   
   progressbar->destroy 
   hist = histogram(ENL_ML,min=0.0,max=20.0,nbins=200)
   hist[0]=0.0
   p = plot(findgen(200)/10.0,hist[0:*])
   ax = p.AXES
   ax[0].TITLE = 'ENL'
   envi_enter_data, ENL_ML, bnames = ['ENL_ML'] 
   print, 'result written to memory, elapsed time: '+strtrim(systime(2)-start,2)+' sec'    
end