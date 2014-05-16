; docformat = 'rst'
; contour_match_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO contour_match_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Contour Matching', $
      REF_VALUE = 'Registration', $
      EVENT_PRO = 'contour_match_run', $
      UVALUE = 'CONTOUR_MATCH', $
      POSITION = 'last'
END

function cc_correlate,contour_a,contour_b
; chain code correlation
   len_a = contour_a.length_cc
   cc_a = contour_a.cc[0:len_a-1]
   len_b = contour_b.length_cc
   cc_b = contour_b.cc[0:len_b-1]
   result=0.0
   for i=0,n_elements(cc_a)-1 do $
      result = correlate(shift(cc_a,i),cc_b) > result
   return,result
end

;+
; :Description:
;       ENVI extension for extraction of tie points for
;       image-image registration.Images may be already
;       georeferenced, in which case GCPs are for "fine
;       adjustement". Uses Laplacian of Gaussian and Sobel
;       filter and contour tracing to match contours::
;          Li, H., Manjunath, B. S., and Mitra, S. K. (1995). 
;          A contour-based approach ;to multisensor image
;          registration. IEEE Transactions on Image Processing,
;          4(3), 320–334.    
; :Params:
;      event:  in, required 
;         if called from the ENVI menu                  
; :Uses:
;       ENVI, CI__DEFINE, COYOTE
; :Author:
;      Mort Canty (2013)       
;-
Pro contour_match_run, event

COMPILE_OPT STRICTARR

corr_thresh = 0.95    ; chain code correlation threshold
length_thresh = 10.0  ; contour length comparison threshold
moment_thresh = 2.0   ; invariant moment comparison threshold
loc_thresh =  100.0   ; position threshold (100 if images georeferenced else 1000000)

   print, '---------------------------------'
   print, 'Contour matching'
   print, systime(0)
   print, '---------------------------------'

; closed contour structure
   closed_contour = { center:fltarr(2), $    ; centroid coordinates
                      moments:fltarr(4),$    ; vector of first 4 Hu moments
                      length_cc:0, $         ; length of chain code
                      best_corr:0.0,$        ; best correlation
                      cc:fltarr(500) }       ; filtered chain code

; GCP table
   n_gcps = 0L
   base_gcps = fltarr(2,2000)
   warp_gcps= fltarr(2,2000)

; get an image band from base image
   envi_select, title='Choose (spatial subset of) base image band',fid=fid_b,dims=dims_b,pos=pos_b,/band_only
   if (fid_b eq -1) then begin
     print,'cancelled'
     return
   endif
   envi_file_query, fid_b, fname=fname_b, xstart=xstart_b, ystart=ystart_b
   print,'Base image chosen: ',fname_b,',   band ',strtrim(pos_b+1,2)
   print,'Dimensions ', dims_b[1:4]
   num_cols = dims_b[2]-dims_b[1]+1
   num_rows = dims_b[4]-dims_b[3]+1

; get an image band from warp image
   envi_select, title='Choose warp image band',fid=fid_w,dims=dims_w,pos=pos_w,/band_only,/no_dims
   if (fid_w eq -1) then begin
     print,'cancelled'
     return
   endif
   envi_file_query, fid_w, fname=fname_w, xstart=xstart_w, ystart=ystart_w

; get LoG sigma
   base = widget_auto_base(title='LoG sigma [*10]')
   wg = widget_sslider(base,  min=10, max=50, $
     value=25, dt=4, uvalue='slide', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
      SIgma = 2.5
      print,'Log Sigma: ',2.5
   end else Sigma = result.slide/10.0


; output file for GCPs
  outfile=dialog_pickfile(filter='*.pts',/write,/overwrite_prompt,title='save CGPs to ASCII')
  if (outfile eq '') then begin
     print,'cancelled'
     return
  endif

; histogram equalize base
   envi_check_save,/transform
   stretch_doit, fid=fid_b,dims=dims_b,pos=[pos_b],method=2,r_fid=r_fid_b,$
     out_min=0,out_max=255, $
     range_by=0,i_min=0.0,i_max=100.0,out_dt=1,/in_memory
   dims1=[-1,0,num_cols-1,0,num_rows-1]
   baseCI = Obj_New('CI',envi_get_data(fid=r_fid_b,dims=dims1,pos=[0]),Sigma)
   envi_file_mng, id=r_fid_b,/remove

; check for georeferencing information
   map_info_b = envi_get_map_info(fid=fid_b)
   map_info_w = envi_get_map_info(fid=fid_w)
   if (map_info_b.mc[2] ne 0.0) and (map_info_w.mc[2] ne 0.0) then begin
; find upper left position of base image within warp image
      envi_convert_file_coordinates,fid_b,dims_b[1],dims_b[3],e,n,/to_map
      envi_convert_file_coordinates,fid_w,X_ul,Y_ul,e,n
      X_ul = round(X_ul)
      Y_ul = round(Y_ul)
; cutout the corresponding spatial subset of the warp image
      dims_w = [-1L,X_ul>dims_w[1],X_ul+num_cols-1<dims_w[2],Y_ul>dims_w[3],Y_ul+num_rows-1<dims_w[4]]
   end else begin
      loc_thresh = 100000.0
      print, 'No georeferencing information available.'
   endelse

   print,'Warp image chosen: ',fname_w,',   band ',strtrim(pos_w+1,2)
   print,'Dimensions ', dims_w[1:4]
   print,'LoG sigma: ',result.slide/10.0

   num_cols = dims_w[2]-dims_w[1]+1
   num_rows = dims_w[4]-dims_w[3]+1

; histogram equalize it
   envi_check_save,/transform
   envi_doit, 'stretch_doit', fid=fid_w,dims=dims_w,pos=[pos_w],method=2,r_fid=r_fid_w,$
     out_min=0,out_max=255, $
     range_by=0,i_min=0.0,i_max=100.0,out_dt=1,/in_memory
   dims1=[-1,0,num_cols-1,0,num_rows-1]
   warpCI = Obj_New('CI',envi_get_data(fid=r_fid_w,dims=dims1,pos=[0]),Sigma)
   envi_file_mng, id=r_fid_w,/remove

  widget_control,/hourglass

; array of raw contours
   max_contours = baseCI->get_max_contours()
   clength = baseCI->get_max_length()
   raw_contours=replicate({ sp:intarr(2),length:0L,closed:-1L,code:bytarr(clength),icode:bytarr(clength) },max_contours)

   print,'Processing base band ...'

   num_raw_contours=0L

   progressbar = Obj_New('cgprogressbar',title='Tracing base contours ...',/cancel)
   progressbar->start
   repeat begin
      if progressbar->CheckCancel() then begin
         print,'cancelled'
         progressbar->Destroy
         return
      endif
      c = baseCI->Trace_Contour()
      if c.closed eq 1 then begin
         raw_contours[num_raw_contours] = c
         num_raw_contours=num_raw_contours+1
         progressbar->Update,((num_raw_contours*100)/max_contours)
      endif
   endrep until (c.closed eq -1) or (num_raw_contours eq max_contours)
   print,num_raw_contours,'    closed contours found in base image'
   progressbar->destroy

;display the base contours (debugging)
;   baseCI->clear_contour_image
;   for i=0,num_raw_contours-1 do baseCI->Write_Contour,raw_contours[i],127
;   envi_enter_data, baseCI->get_contour_image(), xstart=xstart_b+dims_b[1],ystart=ystart_b+dims_b[3]

   widget_control,/hourglass

; generate the set of base contours from the raw contours
   base_contours = replicate(closed_contour,num_raw_contours)
   for i=0,num_raw_contours-1 do begin
     pts = baseCI->to_pixels(raw_contours[i])
     base_contours[i].center[0] = mean(pts[*,0])
     base_contours[i].center[1] = mean(pts[*,1])
     base_contours[i].length_cc = n_elements(pts[*,0])
     base_contours[i].moments = (baseCI->to_moments(raw_contours[i]))[1:4]
     base_contours[i].cc = baseCI->to_filtered_cc(raw_contours[i])
   endfor

; standardize the contour moments
   moments = fltarr(4,num_raw_contours)
   for i=0,num_raw_contours-1 do moments[*,i]=base_contours[i].moments
   for i=0,3 do begin
      sdev = stddev(moments[i,*])
      moments[i,*] = moments[i,*]/sdev
   endfor
   for i= 0,num_raw_contours-1 do base_contours[i].moments = moments[*,i]

   print,'Processing warp band ...'

   num_raw_contours=0L

   progressbar = Obj_New('cgprogressbar',title='Tracing warp contours ...',/cancel)
   progressbar->start
   repeat begin
      if progressbar->CheckCancel() then begin
         print,'cancelled'
         progressbar->Destroy
         return
      endif
      c=warpCI->Trace_Contour()
      if c.closed eq 1 then begin
         raw_contours[num_raw_contours] = c
         num_raw_contours=num_raw_contours+1
         progressbar->Update,((num_raw_contours*100)/max_contours)
      endif
   endrep until (c.closed eq -1) or (num_raw_contours eq max_contours)
   print,num_raw_contours,'    closed contours found in warp image'
   progressbar->destroy

;display the warp contours (debugging)
;   warpCI->clear_contour_image
;   for i=0,num_raw_contours-1 do warpCI->Write_Contour,raw_contours[i],127
;   envi_enter_data, warpCI->get_contour_image(), xstart=xstart_w+dims_w[1],ystart=ystart_w+dims_w[3]

   widget_control,/hourglass

; generate the set of warp contours
   warp_contours = replicate(closed_contour,num_raw_contours)
   for i=0,num_raw_contours-1 do begin
     pts = warpCI->to_pixels(raw_contours[i])
     warp_contours[i].center[0] = mean(pts[*,0])
     warp_contours[i].center[1] = mean(pts[*,1])
     warp_contours[i].length_cc = n_elements(pts[*,0])
     warp_contours[i].moments = (warpCI->to_moments(raw_contours[i]))[1:4]
     warp_contours[i].cc = warpCI->to_filtered_cc(raw_contours[i])
   endfor

; standardize the contour moments
   moments = fltarr(4,num_raw_contours)
   for i=0,num_raw_contours-1 do moments[*,i]=warp_contours[i].moments
   for i=0,3 do begin
      sdev = stddev(moments[i,*])
      moments[i,*] = moments[i,*]/sdev
   endfor
   for i= 0,num_raw_contours-1 do warp_contours[i].moments = moments[*,i]

; ******* contour matching *******

print,'Matching contours ...'

; loop over base_contours
   progressbar = Obj_New('cgprogressbar',title='Matching contours ...',/cancel)
   progressbar->start
   n_base = n_elements(base_contours)
   for i=0L,n_base-1 do begin
      candidates= intarr(n_elements(warp_contours))
      if progressbar->CheckCancel() then begin
         print,'cancelled'
         progressbar->Destroy
         return
      endif
      progressbar->Update,i*100/n_base

; get possible match candidates from moments and center positions
      k=0
      for j=0,n_elements(warp_contours)-1 do begin
         d = abs(base_contours[i].length_cc - warp_contours[j].length_cc)
         dd =  norm(base_contours[i].moments-warp_contours[j].moments)
         ddd = norm(base_contours[i].center-warp_contours[j].center)
         if (d le length_thresh) and (dd le moment_thresh) and (ddd le loc_thresh) then begin
            candidates[k]=j
            k=k+1
         endif
      endfor

; correlate chain codes
      best_corr=0.0
      match = 0
      for j=0,k-1 do begin
         jj = candidates[j]
         corr = cc_correlate(base_contours[i],warp_contours[jj])
         if corr gt best_corr then begin
            best_corr = corr
            match = jj
         endif
      endfor
      if (best_corr gt corr_thresh) and $
         (best_corr gt warp_contours[match].best_corr) then begin
         warp_contours[match].best_corr = best_corr
         base_gcps[*,n_gcps]  = base_contours[i].center
         warp_gcps[*,n_gcps] = warp_contours[match].center
;Debugging
;        warpCI->write_contour,raw_contours[match],255
         n_gcps=n_gcps+1
      endif
   endfor

;Debugging
;   envi_enter_data, warpCI->get_contour_image(),xstart=xstart_w+dims_w[1],ystart=ystart_w+dims_w[3]

   progressbar->Destroy

   print,'Tentative GCPs found: ',strtrim(n_gcps,1)

   print,'Removing outliers ...'

   if n_gcps lt 2 then begin
      print,'No GCPs found'
      void=dialog_message('No GCPs found',/information)
      goto, done
   endif
 ; plausibility check
   ratios = fltarr(n_gcps*(n_gcps-1)/2,3)
   k=0L
   for i=0,n_gcps-1 do for j=i+1,n_gcps-1 do begin
      den = norm(warp_gcps[*,i]-warp_gcps[*,j])
      if den gt 0 then ratios[k,*] = [norm(base_gcps[*,i]-base_gcps[*,j])/den,i,j] $
      else ratios[k,*]=[10,i,j]
      k=k+1
   endfor
   ratios[*,0]=alog(ratios[*,0])

   window,10,xsize=400,ysize=400,title='Ratio cluster'
   wset,10
   hist = histogram(ratios[*,0],nbins=50,min=-1.0,max=1.0,reverse_indices=R)
   plot, hist, color=0, background='FFFFFF'XL
; just take the ratios in the histogram maximum bin
   i = (where(hist eq max(hist)))[0]
   max_indices = R[R[i] : R[i+1]-1]
   ratios = ratios[max_indices,*]
   sigma = 1.0
   k=0L
   n_ratios=n_elements(ratios[*,0])
   while (sigma gt 0.002) and (k lt 5000) do begin
      mn = mean(ratios[*,0])
      indices = where(abs(ratios[*,0]-mn) lt 3*sigma,count)
      if count gt 0 then ratios=ratios[indices,*]
      if n_elements(ratios[*,0]) lt n_ratios then begin
         sigma = stddev(ratios[*,0])
         n_ratios = n_elements(ratios[*,0])
      end else sigma = 2*sigma/3
      k=k+1
   endwhile
   if sigma le 0.002 then begin
      indices = where(abs(ratios[*,0]-mn) lt sigma,count)
      if count gt 0 then ratios=ratios[where(abs(ratios[*,0]-mn) lt sigma),*]
      gcp_table = fltarr(4,n_gcps)
      k=0
;Extract the GCPs from the array of ratios
      for i=0,n_gcps-1 do begin
         indices = where(round(ratios[*,1]) eq i,count1)
         indices = where(round(ratios[*,2]) eq i,count2)
         if (count1 ne 0) and (count2 ne 0) then begin
            gcp_table[0:1,k]=base_gcps[*,i]+[xstart_b,ystart_b]+[dims_b[1],dims_b[3]]
            gcp_table[2:3,k]=warp_gcps[*,i]+[xstart_w,ystart_w]+[dims_w[1],dims_w[3]]
            k=k+1
         endif
      endfor
;Save gcps to 'pts' file
      if k gt 0 then begin
         print, 'GCPs found: ',strtrim(k,2)
         void=dialog_message(strtrim(k,2)+' GCPs found',/information)
         openw,lun,outfile,/get_lun
         printf,lun,';CGPS from contour match'
         printf,lun,';'+systime(0)
         printf,lun,';base file '+fname_b+'  band '+strtrim(pos_b+1,2)
         printf,lun,';Warp file '+fname_w+'  band '+strtrim(pos_w+1,2)
         printf,lun,gcp_table[*,0:k-1]
         free_lun,lun
         print, 'GCPs written to '+outfile
      end else begin
         print,'No GCPs found'
         void=dialog_message('No GCPs found',/information)
      endelse
   end else begin
      print,'No GCPs found'
      void=dialog_message('No GCPs found',/information)
   endelse

   wdelete,10
done:
   print,'done'

   Obj_Destroy, baseCI
   Obj_Destroy, warpCI

End