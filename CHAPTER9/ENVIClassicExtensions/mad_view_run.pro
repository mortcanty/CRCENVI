; docformat = 'rst'
; mad_view_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO mad_view_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'MAD View', $
      REF_VALUE = 'Change Detection Statistics', $
      EVENT_PRO = 'mad_view_run', $
      UVALUE = 'MAD_VIEW',$
      POSITION = 'after'
END

; Utility functions

function Pdist, x, M, V
   return, exp(-(x-M)^2/(2*V))/sqrt(2*!Pi*V)
end

function ProbKS, d
; Kolmogorov-Smirnov probability function (Numerical Recipes, p. 631)
   eps1 = 10^(-6)
   eps2 = 10^(-16)
   fac = 2.0
   sum = 0.0
   termbf = 0.0
   a2 = -2.0*d*d
   for j=1,100 do begin
      term = fac*exp(a2*j*j)
      sum = sum + term
      if (abs(term) le eps1*termbf) or (abs(term) le eps2*sum) then return, sum
      fac=-fac
      termbf = abs(term)
   endfor
   return, 1.0
end

; GUI event handlers

pro mad_thresh_tlb_events
   return
end

pro mad_thresh_cleanup, tlb
   widget_control, tlb, get_Uvalue=info,/no_copy
   if n_elements(info) eq 0 then return
   ptr_free,info.MADs_standard
   ptr_free,info.Chi_Sqr
   ptr_free,info.MADs_byte
   ptr_free,info.map_info
   ptr_free,info.mads_pos
   ptr_free,info.mads_renorm
   print, 'done'
end

pro ONload, event
   widget_control, event.top, get_Uvalue=info,/no_copy
   info.normal_flag=0
   widget_control,info.button_reverseID,set_value='  Normal  '
   envi_select, title='Choose MAD image for thresholding', fid=fid, dims=dims, pos=pos
   if (fid eq -1) then begin
      widget_control,event.top,Set_Uvalue=info,/no_copy
      return
   endif

   widget_control,/hourglass
   info.background=127

   envi_file_query, fid,fname=fname,xstart=xstart, ystart=ystart, $
                    descrip=descrip, bnames=bnames
   info.fname=fname
   info.num_cols = dims[2]-dims[1]+1
   info.num_rows = dims[4]-dims[3]+1
   info.num_pixels = info.num_cols*info.num_rows

; check for chi_sqr band (won't be present if MAD/MNF or if MAD subset was chosen)
   num_bands=n_elements(pos)
   if strmid(bnames[pos[num_bands-1]],0,3) eq 'CHI' then begin
      info.Chi_Sqr = Ptr_New(envi_get_data(fid=fid,dims=dims,pos=pos[num_bands-1]))
      info.num_mads = num_bands-1
   end else begin ; empty array
      info.Chi_Sqr = Ptr_New(fltarr(info.num_cols,info.num_rows))
      info.num_mads = num_bands
   endelse
   if info.num_mads lt 3 then begin
      void=dialog_message('Not enough channels',/error)
      widget_control,event.top,Set_Uvalue=info,/no_copy
      return
   endif
; get IR MAD no-change variances, if present. Not used in calculations.
   vMs = envi_get_header_value(fid, 'varMADs', /float, undefined=undefined)
   if not undefined then $
     info.sigIRMADS = Ptr_New(sqrt(reverse(vMs))) $      ; (envi_get_header_value inverts array order !!!)
   else  info.sigIRMADS = Ptr_New(fltarr(info.num_mads)) ; array of zeros

   info.current_mads[0:2] = pos[0:2]
   info.mads_pos = Ptr_New(pos)
   info.mads_renorm = Ptr_New(intarr(info.num_mads))
   info.map_info = Ptr_New(envi_get_map_info(fid=fid))
   (*info.map_info).mc[0:1]= (*info.map_info).mc[0:1] - [dims[1],dims[3]]
   info.xstart = xstart+dims[1]
   info.ystart = ystart+dims[3]


   if info.num_pixels gt 16000000 then begin
      void=dialog_message('Spatial subset is too large',/error)
      widget_control,event.top,Set_Uvalue=info,/no_copy
      return
   endif
   info.txt=file_basename(info.fname)+',  R('+strtrim(info.current_mads[0]+1,2)+ $
                               ')  G('+strtrim(info.current_mads[1]+1,2)+$
                               ')  B('+strtrim(info.current_mads[2]+1,2)+')'
   widget_control, info.drawID, draw_xsize=info.num_cols, draw_ysize=info.num_rows
   widget_control,info.textID,set_value='Loading ...'
   MADs=fltarr(info.num_mads,info.num_pixels)
   sigMADs = fltarr(info.num_mads)
   for i=0,info.num_mads-1 do begin
      MADs[i,*]=envi_get_data(fid=fid,dims=dims,pos=pos[i])
      mn = mean(MADs[i,*])
      sigMADs[i]=stddev(MADs[i,*])
      MADs[i,*] = (MADs[i,*]-mn)/sigMADs[i]
   endfor

print,sigMADs

   info.sigMADs = Ptr_New(sigMADs)
   info.MADs_standard = Ptr_New(MADs)

   if max(*info.chi_sqr) eq 0 then *info.chi_sqr = reform(total(MADs^2,1),info.num_cols,info.num_rows)

   info.scale=8L
   MADs = ((*info.MADs_standard+8)*255)/16
   MADs = MADs < 255
   MADs = MADs > 0
   info.MADs_byte = Ptr_New(byte(MADs))
   widget_control,info.saveID,sensitive=1
   widget_control,info.assocID,sensitive=1
   widget_control,info.corrID, sensitive=1
   widget_control,info.chi_sqrID, sensitive=1
   widget_control,info.slider_rpID,set_value=32,sensitive=1
   widget_control,info.slider_rmID,set_value=32,sensitive=1
   widget_control,info.slider_gpID,set_value=32,sensitive=1
   widget_control,info.slider_gmID,set_value=32,sensitive=1
   widget_control,info.slider_bpID,set_value=32,sensitive=1
   widget_control,info.slider_bmID,set_value=32,sensitive=1
   widget_control,info.button_ReverseID,set_button=0
   info.thresh_rp = 32
   info.thresh_rm = 32
   info.thresh_gp = 32
   info.thresh_gm = 32
   info.thresh_bp = 32
   info.thresh_bm = 32
   widget_control,info.button_threshID, sensitive=1
   widget_control,info.button_BlackID, sensitive=1
   widget_control,info.button_GrayID, sensitive=1
   widget_control,info.button_ReverseID, sensitive=1
   widget_control,info.button_ChiSqrID, sensitive=1
   widget_control,info.button_RefreshID, sensitive=0
   widget_control,info.slider_scaleID, set_value=8, sensitive=1
   widget_control,info.button_GrayID,/set_button
   widget_control,info.button_BlackID,set_button=0
   widget_control,info.textID,set_value=info.txt
   widget_control,event.top,Set_Uvalue=info,/no_copy
   onButton_dis, event
end

pro ONSlider_scale, event
   widget_control,event.top,Get_Uvalue=info,/no_copy
   oldscale = float(info.scale)
   info.scale = event.value
   newscale = float(info.scale)
   MADs = ((*info.MADs_standard+info.scale)*255)/(2*info.scale)
   MADs = MADs < 255
   MADs = MADs > 0
   *info.MADs_byte = byte(MADs)
   th = 256/(info.scale)
   widget_control,info.slider_rpID,set_value=th,sensitive=1
   widget_control,info.slider_rmID,set_value=th,sensitive=1
   widget_control,info.slider_gpID,set_value=th,sensitive=1
   widget_control,info.slider_gmID,set_value=th,sensitive=1
   widget_control,info.slider_bpID,set_value=th,sensitive=1
   widget_control,info.slider_bmID,set_value=th,sensitive=1
   widget_control,info.button_ReverseID,set_button=0
   info.thresh_rp = th
   info.thresh_rm = th
   info.thresh_gp = th
   info.thresh_gm = th
   info.thresh_bp = th
   info.thresh_bm = th
; reset autothresholds
   cm = info.current_mads
   info.authresh_plus[cm] = info.authresh_plus[cm]*(oldscale/info.scale)
   info.authresh_minus[cm] = info.authresh_minus[cm]*(oldscale/info.scale)
   widget_control,info.button_RefreshID, sensitive=1
   widget_control,event.top,set_Uvalue=info,/no_copy
end

pro ONquit, event
   widget_control,event.top,get_Uvalue=info,/no_copy
   ptr_free,info.MADs_byte
   ptr_free,info.MADs_standard
   ptr_free,info.map_info
   ptr_free,info.mads_pos
   widget_control, event.top, /destroy
   print,'done'
end

pro ONbutton_thresh, event
   widget_control,event.top,Get_Uvalue=info,/no_copy
   scale = info.scale
; loop over MAD components
   widget_control,info.textID,set_value='Running EM algorithm ...'
   widget_control,/hourglass
   indices = randomu(seed,10000,/long) mod info.num_pixels
   for i=0,2 do begin
      k = (where(*info.MADs_pos eq info.current_mads[i]))[0]
      cm = info.current_mads[i]
      if (*info.mads_renorm)[k] eq 0 then begin
;         S=transpose((*info.MADs_standard)[k,indices])
         S=(*info.MADs_standard)[k,indices]
; run the EM algorithm
         unfrozen = where((abs(S) gt 0.02) and (S lt 2.0) and (S gt -2.0),complement=frozen)
         unfrozen = where( (S lt 2.0) and (S gt -2.0),complement=frozen)
         U = randomu(seed,n_elements(S),3)
         U[frozen,*] = 0.0
         ind = where(S le -2.0,count)
         if count gt 0 then U[ind,0] = 1.0
         ind = where(abs(S) le 0.02,count)
         if count gt 0 then U[ind] = 1.0
         ind = where(S ge 2.0,count)
         if count gt 0 then U[ind,2] = 1.0
         for j=0,2 do U[*,j]=U[*,j]/total(U,2)
         em, S, U, Ms, Ps, Vs, unfrozen=unfrozen, T0=0.0
         Mcm = Ms[0]
         Mnc = Ms[1]
         Mcp = Ms[2]
         Vcm = Vs[0]
         Vnc = Vs[1]
         Vcp = Vs[2]
         Pcm = Ps[0]
         Pnc = Ps[1]
         Pcp = Ps[2]
; calculate the change thresholds and correct to the new scale, i.e. multiply by (256/(2*scale))/sqrt(Vnc), see below
         A  = alog (sqrt(Vnc/Vcp)*(Pcp/Pnc))
         info.authresh_plus[cm] = (Mcp*Vnc-Mnc*Vcp-sqrt(Vnc*Vcp)*sqrt((Mnc-Mcp)^2+2*A*(Vnc-Vcp)))/(Vnc-Vcp)*(256/(2*scale))/sqrt(Vnc)
         A  = alog (sqrt(Vnc/Vcm)*(Pcm/Pnc))
         info.authresh_minus[cm] = (Mcm*Vnc-Mnc*Vcm+sqrt(Vnc*Vcm)*sqrt((Mnc-Mcm)^2+2*A*(Vnc-Vcm)))/(Vnc-Vcm)*(256/(2*scale))/sqrt(Vnc)
; plot fits
         if info.normal_flag eq 0 then h = histogram((*info.MADs_byte)[k,*],min=0,max=255) $
                             else h = histogram(255-(*info.MADs_byte)[k,*],min=0,max=255)
         x= (findgen(256)-127.0)/(256/(2*scale))
         pcm1 = Pcm*pdist(x,Mcm,Vcm)/(256/(2*scale))
         pnc1 = Pnc*pdist(x,Mnc,Vnc)/(256/(2*scale))
         pcp1 = Pcp*pdist(x,Mcp,Vcp)/(256/(2*scale))
         
;         pp = plot(findgen(256),alog10([h]),thick= 2,font_name="times", font_size=20)
;         pp = plot(findgen(256),alog10([total(h)*pcm1]+1),/overplot,thick= 2)
;         pp = plot(findgen(256),alog10([total(h)*pnc1]+1),/overplot,thick= 2)
;         pp = plot(findgen(256),alog10([total(h)*pcp1]+1),/overplot,thick= 2)
;         pp = plot(findgen(256),alog10([total(h)*(pnc1+pcm1+pcp1)]+1),/overplot)
         
         envi_plot_data,findgen(256),alog10([[h],$
                                   [total(h)*pcm1],$
                                   [total(h)*pnc1],$
                                   [total(h)*pcp1],$
                                   [total(h)*(pnc1+pcm1+pcp1)]]+1),$
                                   plot_styles=[1,0,0,0,0],$
                                   title='MAD '+strtrim(info.current_mads[i]+1,2),$
                                   ytitle='log(counts)',$
                                   base=base
         widget_control, base, group_leader=event.top
; recalculate the standardized distribution for the current MAD components
; using the no-change variances determined from the EM algorithm
; ie: s1^2 = first estimate of no-change variance, MADs_standard = MADs_raw/s1
;     s2^2 = better estimate of no-change variance after threshold calculation
;     corrected MADs_standard = MADs_raw/s2 = MADs_standard*s1/s2 = MADs_standard/(s2/s1)
;     and s2/s1 = sqrt(Vnc)
         s1 = (*info.sigMADs)[k]
         (*info.MADs_standard)[k,*] = (*info.MADs_standard)[k,*]/sqrt(Vnc)
         print,'MAD '+strtrim((*info.mads_pos)[k]+1,2)
         print,'  sig(MAD)            : '+strtrim(s1,2)
         print,'  sig(MAD) (no change): '+strtrim(s1*sqrt(Vnc),2)
         print,'  sqrt(2(1-rho))      : '+strtrim((*info.sigIRMADs)[k],2)
         print,'  No change fraction  : '+strtrim(Pnc,2)
         (*info.mads_renorm)[k] = 1
      endif
      if i eq 0 then begin
         info.thresh_rp = fix(info.authresh_plus[cm])
         info.thresh_rm = -fix(info.authresh_minus[cm])
         widget_control,info.slider_rpID,set_value=info.thresh_rp
         widget_control,info.slider_rmID,set_value=info.thresh_rm
      end else if i eq 1 then begin
         info.thresh_gp = fix(info.authresh_plus[cm])
         info.thresh_gm = -fix(info.authresh_minus[cm])
         widget_control,info.slider_gpID,set_value=info.thresh_gp
         widget_control,info.slider_gmID,set_value=info.thresh_gm
      end else if i eq 2 then begin
         info.thresh_bp = fix(info.authresh_plus[cm])
         info.thresh_bm = -fix(info.authresh_minus[cm])
         widget_control,info.slider_bpID,set_value=info.thresh_bp
         widget_control,info.slider_bmID,set_value=info.thresh_bm
      end
   endfor
; re-stretch the histograms
   MADs = ((*info.MADs_standard+scale)*255)/(2*scale)
   MADs = MADs < 255
   MADs = MADs > 0
   *info.MADs_byte = byte(MADs)
   if info.normal_flag then *info.MADs_byte=255-*info.MADs_byte
   widget_control,info.button_RefreshID, sensitive=1
   widget_control,info.textID,set_value=info.txt
   widget_control,event.top,Set_Uvalue=info,/no_copy
end

pro ONassoc, event
   widget_control,event.top,get_Uvalue=info,/no_copy
   txt = ' '
   for i=0,info.num_mads-1 do txt=txt+strtrim((*info.mads_pos)[i]+1,2)+' '
   base = widget_auto_base(title='Available MADs:'+txt)
   we = widget_edit(base,  list=['R','G','B'], uvalue='edit', $
        vals=info.current_mads+1, dt=2, /auto)
   result = auto_wid_mng(base)
   if (result.accept ne 0) then begin
      if (where(result.edit[0]-1 eq *info.MADs_pos) eq -1) or $
         (where(result.edit[1]-1 eq *info.MADs_pos) eq -1) or $
         (where(result.edit[2]-1 eq *info.MADs_pos) eq -1) then begin
         void = dialog_message('Invalid MAD index',/error)
         widget_control,event.top,set_Uvalue=info,/no_copy
         return
      endif
      current_mads = result.edit-1
      info.txt=file_basename(info.fname)+',  R('+strtrim(current_mads[0]+1,2)+ $
                               ')  G('+strtrim(current_mads[1]+1,2)+$
                               ')  B('+strtrim(current_mads[2]+1,2)+')'
      print,info.txt
      info.current_mads = current_mads
      widget_control,info.textID,set_value=info.txt
      widget_control,info.slider_rpID,set_value=256/info.scale,sensitive=1
      widget_control,info.slider_rmID,set_value=256/info.scale,sensitive=1
      widget_control,info.slider_gpID,set_value=256/info.scale,sensitive=1
      widget_control,info.slider_gmID,set_value=256/info.scale,sensitive=1
      widget_control,info.slider_bpID,set_value=256/info.scale,sensitive=1
      widget_control,info.slider_bmID,set_value=256/info.scale,sensitive=1
      info.thresh_rp = 256/info.scale
      info.thresh_rm = 256/info.scale
      info.thresh_gp = 256/info.scale
      info.thresh_gm = 256/info.scale
      info.thresh_bp = 256/info.scale
      info.thresh_bm = 256/info.scale
      widget_control,info.button_RefreshID, sensitive=1
   endif
   widget_control,event.top,set_Uvalue=info,/no_copy
end

pro ONcorr, event
   widget_control,event.top,get_Uvalue=info,/no_copy
; mask the change pixels
   MADs = (*info.MADs_standard)[info.current_mads,*]

   cmask= where( ((*info.MADs_Byte)[0,*] gt 127+info.thresh_rp) or $
                 ((*info.MADs_Byte)[0,*] lt 127-info.thresh_rm) or $
                 ((*info.MADs_Byte)[1,*] gt 127+info.thresh_gp) or $
                 ((*info.MADs_Byte)[1,*] lt 127-info.thresh_gm) or $
                 ((*info.MADs_Byte)[2,*] gt 127+info.thresh_bp) or $
                 ((*info.MADs_Byte)[2,*] lt 127-info.thresh_bm), count)

; get image to be correlated
   envi_select, title='Choose image for correlation with MAD', fid=fid, dims=dims,pos=pos
   if fid eq -1 then begin
      widget_control,event.top,Set_Uvalue=info,/no_copy
      return
   endif
   if (dims[2]-dims[1]+1 ne info.num_cols) or (dims[4]-dims[3]+1 ne info.num_rows) then begin
      m=dialog_message('Dimensions dont match',/error)
      widget_control,event.top,Set_Uvalue=info,/no_copy
      return
   endif
   envi_file_query, fid,fname=fname
   num_bands = n_elements(pos)
   mads_change  = fltarr(3,n_elements(cmask))
   image_change = fltarr(num_bands,n_elements(cmask))
   for i=0,2 do mads_change[i,*]=MADs[i,cmask]
   for i=0,num_bands-1 do begin
      band = envi_get_data(fid=fid,dims=dims,pos=pos[i])
      image_change[i,*] = band[cmask]
   endfor
; correlation
   samples = transpose( [[transpose(mads_change)],[transpose(image_change)]] )
   print, 'correlation of MADs ',info.current_mads+1,' with original bands from '+fname
   print, n_elements(cmask),' change pixels only'
   corr = (correlate(samples,/double))[3:num_bands+2,0:2]
   print,corr
; plot results (first row is white, then rgb)
   arr = [[fltarr(num_bands)],[corr]]
   
;   ppp = plot( findgen(num_bands)+1, arr[*,1],thick= 2,font_name="times", font_size=26)

   
   envi_plot_data, findgen(num_bands)+1, arr, title='Correlations with '+file_basename(fname), xtitle='band', $
     plot_names = [' ','MAD'+strtrim(info.current_mads[0]+1,1),'MAD'+strtrim(info.current_mads[1]+1,1), $
                   'MAD'+strtrim(info.current_mads[2]+1,1)], $
     plot_styles = [0,0,2,3], $
     base=base
   widget_control, base, group_leader=event.top
   widget_control,event.top,set_Uvalue=info,/no_copy
end

pro ONchi_sqr, event
   widget_control,event.top,Get_Uvalue=info,/no_copy
   p999=chisqr_cvf(0.001,info.num_mads)
; true chi-sqr distribution, normalized to num_pixels
   X = randomu(seed,info.num_mads,info.num_pixels,/normal)
   hst1 = histogram(total(X^2,1), min=0.0, max=p999, nbins=100)
; chi_sqr at read in time
   hst2 = histogram(*info.Chi_Sqr,  min=0.0, max=p999, nbins=100)
; chi_sqr using standardized MADs
   chi_sqr = total((*info.MADs_standard)^2,1)
   hst3 = histogram(chi_sqr, min=0.0, max=p999, nbins=100)

   envi_plot_data,findgen(100)*p999/100,[[hst1],[hst2],[hst3]], $
                                   plot_styles=[0,0,0], $
                                   plot_names= ['chisqr,'+strtrim(info.num_mads,2),'on load', 'current'], $
                                   title='Chi-sqr distribution',$
                                   base=base
   widget_control, base, group_leader=event.top
   widget_control,event.top,set_Uvalue=info,/no_copy
end

pro ONbutton_Black, event
   widget_control,event.top,Get_Uvalue=info,/no_copy
   widget_control,/hourglass
   info.Background = 0
   widget_control,event.top,set_Uvalue=info,/no_copy
   onbutton_dis,event
end

pro ONbutton_Gray, event
   widget_control,event.top,Get_Uvalue=info,/no_copy
   widget_control,/hourglass
   info.Background = 127
   widget_control,event.top,set_Uvalue=info,/no_copy
   onbutton_dis,event
end

pro ONbutton_Reverse, event
   widget_control,event.top,Get_Uvalue=info,/no_copy
   widget_control,/hourglass
   if info.normal_flag eq 0 then begin
      info.normal_flag=1
      widget_control,info.button_reverseID,set_value='Reversed'
   end else begin
      info.normal_flag=0
      widget_control,info.button_reverseID,set_value='  Normal  '
   endelse
   *info.MADs_byte = 255-*info.MADs_byte
   widget_control,event.top,set_Uvalue=info,/no_copy
   onbutton_dis,event
end

pro ONbutton_ChiSqr, event
   widget_control,event.top,Get_Uvalue=info,/no_copy
   info.chiSqr_flag=1
   wset,info.wID
   chi_sqr = reform(total((*info.MADs_standard)^2,1),info.num_cols,info.num_rows)
   tvscl, chi_sqr
   widget_control,info.button_RefreshID, sensitive=1
   widget_control,event.top,set_Uvalue=info,/no_copy
end

pro ONbutton_Dis, event
   widget_control,event.top,Get_Uvalue=info,/no_copy
   widget_control,/hourglass
   i0 = where(*info.MADs_pos eq info.current_mads[0])
   i1 = where(*info.MADs_pos eq info.current_mads[1])
   i2 = where(*info.MADs_pos eq info.current_mads[2])
   MADs=(*info.MADs_byte)[[i0,i1,i2],*]
   if info.normal_flag eq 0 then begin
      thresh_rp=info.thresh_rp
      thresh_rm=info.thresh_rm
      thresh_gp=info.thresh_gp
      thresh_gm=info.thresh_gm
      thresh_bp=info.thresh_bp
      thresh_bm=info.thresh_bm
   end else begin
      thresh_rp=info.thresh_rm
      thresh_rm=info.thresh_rp
      thresh_gp=info.thresh_gm
      thresh_gm=info.thresh_gp
      thresh_bp=info.thresh_bm
      thresh_bm=info.thresh_bp
   endelse
   nc_mask = where( (MADs[0,*]-127 lt thresh_rp) and (MADs[0,*]-127 gt -thresh_rm), count )
   if count ne 0 then MADs[0,nc_mask]=info.Background
   if thresh_rp eq 128 then begin
      indices = where(MADs[0,*]-127 gt 0, count)
      if count ne 0 then MADs[0,indices]=info.Background
   endif
   if thresh_rm eq 128 then begin
      indices = where(MADs[0,*]-127 lt 0, count)
      if count ne 0 then MADs[0,indices]=info.Background
   endif
   nc_mask = where( (MADs[1,*]-127 lt thresh_gp) and (MADs[1,*]-127 gt -thresh_gm), count )
   if count ne 0 then MADs[1,nc_mask]=info.Background
   if thresh_gp eq 128 then begin
      indices = where(MADs[1,*]-127 gt 0, count)
      if count ne 0 then MADs[1,indices]=info.Background
   endif
   if thresh_gm eq 128 then begin
      indices = where(MADs[1,*]-127 lt 0, count)
      if count ne 0 then MADs[1,indices]=info.Background
   endif
   nc_mask = where( (MADs[2,*]-127 lt thresh_bp) and (MADs[2,*]-127 gt -thresh_bm), count )
   if count ne 0 then MADs[2,nc_mask]=info.Background
   if thresh_bp eq 128 then begin
      indices = where(MADs[2,*]-127 gt 0, count)
      if count ne 0 then MADs[2,indices]=info.Background
   endif
   if thresh_bm eq 128 then begin
      indices = where(MADs[2,*]-127 lt 0, count)
      if count ne 0 then MADs[2,indices]=info.Background
   endif
   wset,info.wID
   tv,reform(MADs,3,info.num_cols,info.num_rows),/true
   widget_control,info.button_RefreshID, sensitive=0
   widget_control,event.top,set_Uvalue=info,/no_copy
end

pro ONslider_rp, event
   widget_control,event.top,get_Uvalue=info,/no_copy
   info.thresh_rp = event.value
   if info.current_mads[0] eq info.current_mads[1] then begin
       info.thresh_gp = event.value
       widget_control,info.slider_gpID,set_value=event.value
   endif
   if info.current_mads[0] eq info.current_mads[2] then begin
       info.thresh_bp = event.value
       widget_control,info.slider_bpID,set_value=event.value
   endif
   widget_control,info.button_RefreshID, sensitive=1
   widget_control,event.top,set_Uvalue=info,/no_copy
end

pro ONslider_rm, event
   widget_control,event.top,get_Uvalue=info,/no_copy
   info.thresh_rm = event.value
   if info.current_mads[0] eq info.current_mads[1] then begin
       info.thresh_gm = event.value
       widget_control,info.slider_gmID,set_value=event.value
   endif
   if info.current_mads[0] eq info.current_mads[2] then begin
       info.thresh_bm = event.value
       widget_control,info.slider_bmID,set_value=event.value
   endif
   widget_control,info.button_RefreshID, sensitive=1
   widget_control,event.top,set_Uvalue=info,/no_copy
end


pro ONslider_gp, event
   widget_control,event.top,get_Uvalue=info,/no_copy
   info.thresh_gp = event.value
   if info.current_mads[1] eq info.current_mads[0] then begin
       info.thresh_rp = event.value
       widget_control,info.slider_rpID,set_value=event.value
   endif
   if info.current_mads[1] eq info.current_mads[2] then begin
       info.thresh_bp = event.value
       widget_control,info.slider_bpID,set_value=event.value
   endif
   widget_control,info.button_RefreshID, sensitive=1
   widget_control,event.top,set_Uvalue=info,/no_copy
end


pro ONslider_gm, event
   widget_control,event.top,get_Uvalue=info,/no_copy
   info.thresh_gm = event.value
   if info.current_mads[1] eq info.current_mads[0] then begin
       info.thresh_rm = event.value
       widget_control,info.slider_rmID,set_value=event.value
   endif
   if info.current_mads[1] eq info.current_mads[2] then begin
       info.thresh_bm = event.value
       widget_control,info.slider_bmID,set_value=event.value
   endif
   widget_control,info.button_RefreshID, sensitive=1
   widget_control,event.top,set_Uvalue=info,/no_copy
end


pro ONslider_bp, event
   widget_control,event.top,get_Uvalue=info,/no_copy
   info.thresh_bp = event.value
   if info.current_mads[2] eq info.current_mads[0] then begin
       info.thresh_rp = event.value
       widget_control,info.slider_rpID,set_value=event.value
   endif
   if info.current_mads[2] eq info.current_mads[1] then begin
       info.thresh_gp = event.value
       widget_control,info.slider_gpID,set_value=event.value
   endif
   widget_control,info.button_RefreshID, sensitive=1
   widget_control,event.top,set_Uvalue=info,/no_copy
end


pro ONslider_bm, event
   widget_control,event.top,get_Uvalue=info,/no_copy
   info.thresh_bm = event.value
   if info.current_mads[2] eq info.current_mads[0] then begin
       info.thresh_rm = event.value
       widget_control,info.slider_rmID,set_value=event.value
   endif
   if info.current_mads[2] eq info.current_mads[1] then begin
       info.thresh_gm = event.value
       widget_control,info.slider_gmID,set_value=event.value
   endif
   widget_control,info.button_RefreshID, sensitive=1
   widget_control,event.top,set_Uvalue=info,/no_copy
end


pro ONsave, event
   widget_control,event.top,get_Uvalue=info,/no_copy
   if info.chiSqr_flag then begin
      info.chiSqr_flag = 0
      envi_enter_data, *info.Chi_sqr + 0, $
                      map_info=*(info.map_info), $
                      xstart=info.xstart, ystart=info.ystart
   end else begin
      MADs=(*info.MADs_byte)[info.current_mads,*]
      nc_mask = where( (MADs[0,*]-127 lt info.thresh_rp) and (MADs[0,*]-127 gt -info.thresh_rm), count )
      if count ne 0 then MADs[0,nc_mask]=info.Background
      if info.thresh_rp eq 128 then begin
         indices = where(MADs[0,*]-127 gt 0, count)
         if count ne 0 then MADs[0,indices]=info.Background
      endif
      if info.thresh_rm eq 128 then begin
         indices = where(MADs[0,*]-127 lt 0, count)
         if count ne 0 then MADs[0,indices]=info.Background
      endif
      nc_mask = where( (MADs[1,*]-127 lt info.thresh_gp) and (MADs[1,*]-127 gt -info.thresh_gm), count )
      if count ne 0 then MADs[1,nc_mask]=info.Background
      if info.thresh_gp eq 128 then begin
         indices = where(MADs[1,*]-127 gt 0, count)
         if count ne 0 then MADs[1,indices]=info.Background
      endif
      if info.thresh_gm eq 128 then begin
         indices = where(MADs[1,*]-127 lt 0, count)
         if count ne 0 then MADs[1,indices]=info.Background
      endif
      nc_mask = where( (MADs[2,*]-127 lt info.thresh_bp) and (MADs[2,*]-127 gt -info.thresh_bm), count )
      if count ne 0 then MADs[2,nc_mask]=info.Background
      if info.thresh_bp eq 128 then begin
         indices = where(MADs[2,*]-127 gt 0, count)
         if count ne 0 then MADs[2,indices]=info.Background
      endif
      if info.thresh_bm eq 128 then begin
         indices = where(MADs[2,*]-127 lt 0, count)
         if count ne 0 then MADs[2,indices]=info.Background
      endif
      MADs1 = bytarr(info.num_cols,info.num_rows,3)
      for i=0,2 do MADs1[*,*,i]=reform(MADs[i,*],info.num_cols,info.num_rows)
      envi_enter_data,MADs1, $
                      map_info=*(info.map_info), $
                      xstart=info.xstart, ystart=info.ystart
   endelse
   widget_control,event.top,set_Uvalue=info,/no_copy
end

pro onHelp, event
   help = file_which('madViewHelp.pdf',/include_current_dir)
   ONLINE_HELP, 'MadView', BOOK=help
end

PRO onAbout, event
   dummy=dialog_message(['MAD_VIEW Version 1.0 ',$
                         'M. Canty, 2013'],/info)
END

;+
; :Description:
;       GUI for viewing and thresholding MAD images  
; :Params:
;       event:  in, optional 
;          required if called from ENVI              
; :Uses:
;       ENVI::       
;       EM::       
;       COYOTE
; :Author:
;       Mort Canty (2013)        
;-
PRO mad_view_run, event

print, '-------------------------'
print, 'MAD View'
print, systime(0)
print, '-------------------------'

; widget creation and registration

tlb           = widget_base(/column,title="MAD View",mbar=menubarID,xoffset=500,yoffset=100)
fileID        = widget_button(menubarID,value='File',/menu)
dataID        = widget_button(menubarID,value='Data',/menu)
helpID        = widget_button(menubarID,value='Help',/menu)
mvHelpID      = widget_button(helpID,value='MAD View Help',event_pro='ONhelp')
mvAboutID     = widget_button(helpID,value='About MAD View',event_pro='ONabout')
loadID        = widget_button (fileID, value='Load', event_pro='ONload')
saveID        = widget_button(fileID,value='To ENVI',event_pro='ONsave',Sensitive=0)
quitID        = widget_button(fileID,value='Quit',event_pro='ONquit')
assocID       = widget_button(dataID,value='Associate',event_pro='ONassoc',Sensitive=0)
corrID        = widget_button(dataID,value='Correlate',event_pro='ONcorr',Sensitive=0)
chi_sqrID     = widget_button(dataID,value='Chi-Sqr',event_pro='ONchi_sqr',Sensitive=0)
textID        = widget_text(tlb,value=' ')
slider_top    = widget_base(tlb,/row,ysize=70)
slider_top_rm = widget_base(slider_top,/column,ysize=70)
slider_rmID   = widget_slider(slider_top_rm,xsize=100,value=32,min=0,max=128,event_pro='ONslider_rm',Sensitive=0)
void          = widget_text(slider_top_rm,value='   red negative')
slider_top_rp = widget_base(slider_top,/column)
slider_rpID   = widget_slider(slider_top_rp,xsize=100,value=32,min=0,max=128,event_pro='ONslider_rp',Sensitive=0)
void          = widget_text(slider_top_rp,value='   red positive')
slider_top_gm = widget_base(slider_top,/column)
slider_gmID   = widget_slider(slider_top_gm,xsize=100,value=32,min=0,max=128,event_pro='ONslider_gm',Sensitive=0)
void          = widget_text(slider_top_gm,value='   green negative')
slider_top_gp = widget_base(slider_top,/column)
slider_gpID   = widget_slider(slider_top_gp,xsize=100,value=32,min=0,max=128,event_pro='ONslider_gp',Sensitive=0)
void          = widget_text(slider_top_gp,value='   green positive')
slider_top_bm = widget_base(slider_top,/column)
slider_bmID   = widget_slider(slider_top_bm,xsize=100,value=32,min=0,max=128,event_pro='ONslider_bm',Sensitive=0)
void          = widget_text(slider_top_bm,value='   blue negative')
slider_top_bp = widget_base(slider_top,/column)
slider_bpID   = widget_slider(slider_top_bp,xsize=100,value=32,min=0,max=128,event_pro='ONslider_bp',Sensitive=0)
void          = widget_text(slider_top_bp,value='   blue positive')

button_top = widget_base(tlb,/row,ysize=40)
button_ReverseID = widget_button(button_top, value= '  Normal  ',event_pro='Onbutton_Reverse',Sensitive=0)
button_threshID  = widget_button(button_top, value='  Auto thresh  ',event_pro='ONbutton_thresh',Sensitive=0)
button_RefreshID  = widget_button(button_top, value= '     Refresh     ',event_pro='Onbutton_Dis',Sensitive=0)
button_ChiSqrID = widget_button(button_top, value= '     Chi Square     ',event_pro='Onbutton_ChiSqr',Sensitive=0)
slider_top_scale = widget_base(button_top,/column)
slider_scaleID   = widget_slider(slider_top_scale,xsize=100,value=8,min=4,max=16,event_pro='ONslider_scale',Sensitive=0)
radio_top   = widget_base(button_top,/row,/exclusive)
button_BlackID= widget_button(radio_top, value='Black     ',event_pro='Onbutton_Black',Sensitive=0)
button_GrayID = widget_button(radio_top, value= 'Gray       ',event_pro='Onbutton_Gray',Sensitive=0)
drawID        = widget_draw(tlb,xsize=800,ysize=700, X_Scroll_size=800, Y_Scroll_size=700,/scroll)
widget_control, tlb, /realize
widget_control, drawID, get_value=wID

; widget communication structure

info =  {  MADs_byte: Ptr_New(), $     ; MAD (or MAD/MNF) variates in byte streched initially to +- 8 stdvs
           MADs_standard:Ptr_New(), $  ; standardized MAD variates in float
           Chi_Sqr: Ptr_New(), $       ; Chi-square image from MAD calculation
           sigMADs: Ptr_New(), $       ; standard deviations of no-change pixels from 3-Gaussian fits
           sigIRMADs: Ptr_New(), $     ; standard deviations of no-change IR-MADs (2(1-rho)), set by MAD_RUN)
           map_info:  Ptr_New(), $     ; map info
           mads_pos:  Ptr_New(), $     ; MAD variates read in
           mads_renorm: Ptr_New(), $   ; renormalized MAD variates
           authresh_minus: fltarr(50), $
           authresh_plus: fltarr(50), $
           fname: ' ',$
           txt: ' ',$
           xstart: 1L,$
           ystart: 1L,$
           textID: textID, $
           saveID: saveID, $
           assocID: assocID, $
           corrID: corrID,$
           chi_sqrID: chi_sqrID, $
           slider_rpID: slider_rpID, $
           slider_rmID: slider_rmID, $
           slider_gpID: slider_gpID, $
           slider_gmID: slider_gmID, $
           slider_bpID: slider_bpID, $
           slider_bmID: slider_bmID, $
           drawID: drawID, $
           wID:      wID, $
           button_threshID:button_threshID,$
           button_BlackID:button_BlackID,$
           button_GrayID:button_GrayID,$
           button_ReverseID:button_ReverseID,$
           button_RefreshID:button_RefreshID,$
           button_ChiSqrID:button_ChiSqrID,$
           slider_scaleID:slider_scaleID,$
           current_mads: [0,1,2], $
           num_cols: 500L, $
           num_rows: 500L, $
           num_pixels: 250000L, $
           num_mads: 3,   $
           chiSqr_flag: 0, $
           scale: 16L,    $
           thresh_rp: 16L, $
           thresh_rm: 16L, $
           thresh_gp: 16L, $
           thresh_gm: 16L, $
           thresh_bp: 16L, $
           thresh_bm: 16L, $
           normal_flag: 0, $
           background: 127 $
       }

!ORDER = 1

widget_control,tlb,set_Uvalue=info,/no_copy

XManager, 'mad_thresh',tlb,/no_block, $
           event_handler='mad_thresh_tlb_events', $
           cleanup='mad_thresh_cleanup'

END
