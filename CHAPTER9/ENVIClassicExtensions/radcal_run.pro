; docformat = 'rst'
; radcal_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO radcal_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Radiometric Normalization', $
      REF_VALUE = 'Change Detection Statistics', $
      EVENT_PRO = 'radcal_run', $
      UVALUE = 'RADCAL',$
      POSITION = 'after', $
      /SEPARATOR
END

;+
; :Description:
;     Radiometric calibration using IR-MAD::
;        Canty, M. J. and Nielsen, A. A. (2008). 
;        Automatic radiometric normalization of multitemporal
;        satellite imagery with the iteratively re-weighted MAD
;        transformation. Remote Sensing of Environment, 
;        112(3), 1025–1036.
; :Params:
;       event:  in, optional 
;          required if called from ENVI        
; :Uses:
;       ENVI::       
;       ORTHO_REGRESS::             
;       COYOTE         
; :Author:
;       Mort Canty (2009)  
;-
pro radcal_run, event

COMPILE_OPT STRICTARR

print, '---------------------------------'
print, 'Radiometric calibration with MAD'
print, systime(0)


catch, theError
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
   return
endif

envi_select, title='Choose (subset of) reference image', fid=fid1, dims=dims1,pos=pos1
    if (fid1 eq -1) then begin
      print,'cancelled'
      return
    endif
envi_file_query, fid1, fname=fname1, ns=ns, nl=nl, nb=nb
print, 'Reference: ',fname1

envi_select, title='Choose (subset of) target image', fid=fid2, dims=dims2,pos=pos2
    if (fid2 eq -1) then begin
       print,'cancelled'
       return
    endif
envi_file_query, fid2, fname=fname2, bnames=bnames2, ns=ns, nl=nl, nb=nb, xstart=xstart2, ystart=ystart2, wl=wl2
print, 'Target:    ',fname2
print, '---------------------------------'
bnames2 = 'renorm('+bnames2[pos2]+')'
wl2 = wl2[pos2]

if ((dims1[2]-dims1[1]) ne (dims2[2]-dims2[1])) or $
   ((dims1[4]-dims1[3]) ne (dims2[4]-dims2[3])) or $
   (n_elements(pos1) ne n_elements(pos2)) $
then begin
   print, 'Dimensions of the two images are different. Aborting.'
   Message, 'Dimensions are different. Aborting.'
   return
endif

envi_select, title='Choose (spatial subset of) chi_square image', $
 fid=fid3, dims=dims3,pos=pos3,/band_only,/mask,m_fid=m_fid
    if (fid3 eq -1) then begin
      print,'cancelled'
      return
    endif
envi_file_query, fid3, fname=fname3, ns=ns, nl=nl, nb=nb3
print, 'chi-square image: ',fname3

if ((dims1[2]-dims1[1]) ne (dims3[2]-dims3[1])) or $
   ((dims1[4]-dims1[3]) ne (dims3[4]-dims3[3])) then begin
   print, 'Dimensions of the two images are different. Aborting.'
   Message, 'Dimensions are different. Aborting.'
endif

num_cols = dims1[2]-dims1[1]+1
num_rows = dims1[4]-dims1[3]+1
num_pixels = (num_cols*num_rows)
num_bands = n_elements(pos1)
if m_fid ne -1 then mask = envi_get_data(fid=m_fid,dims=dims1,pos=0) $
      else mask = bytarr(num_cols,num_rows)+1B
      
; read and plot chi-square histogram out to 3 times 99.9 percentile
p999=chisqr_cvf(0.001,nb3-1)
chi_sqr = envi_get_data(fid=fid3,dims=dims3,pos=pos3)
ncp = 1.0 - chisqr_pdf(chi_sqr,nb3-1)
idx = where(mask,count,complement=idxc,ncomplement=ncomplement)
if ncomplement gt 0 then ncp[idxc] = 0.0
X = randomu(seed,nb3-1,n_elements(idx),/normal)
X2 = total(X^2,1)
hst = histogram(X2, min=0.0, max=3*p999, nbins=500)
hst1 = histogram(chi_sqr[idx], min=0.0, max=3*p999, nbins=500)
envi_plot_data,findgen(500)*3*p999/500,[[hst],[hst1]], $
                                   title='Chi-Sqr Histogram',$
                                   ytitle='Frequency', $
                                   xtitle='chi_sqr', $
                                   xoff=100, yoff=100

; select no-change
base = widget_auto_base(title='No-change threshold')
wg = widget_sslider(base, title='Select minimum no-change probability (% x 10)', min=50, max=999, scale=1, $
   value=950, dt=1, uvalue='slide', xs=[300,7],/auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print,'cancelled'
   return
endif
prct = float(result.slide)/1000.0
mask_nochange = where(ncp ge prct)


; output destination
base = widget_auto_base(title='Save normalized image')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result = auto_wid_mng(base)

widget_control,/hourglass

; tie point
map_info2 = envi_get_map_info(fid=fid2)
envi_convert_file_coordinates, fid2, dims2[1], dims2[3], e, n, /to_map
map_info2.mc= [0D,0D,e,n]

image3=fltarr(num_bands,num_pixels)

mrs=fltarr(num_bands) ; means of reference image
mts=fltarr(num_bands) ; means of target image
mns=fltarr(num_bands) ; means of normalized image
vrs=fltarr(num_bands) ; variances of reference image
vts=fltarr(num_bands) ; variances of target image
vns=fltarr(num_bands) ; variances of normalized image
df=intarr(num_bands) ; degrees of freedom
aa= fltarr(num_bands) ; slope of orthogonal regression curve
xm = fltarr(num_bands); mean of X
ym = fltarr(num_bands); mean of Y
table = fltarr(9,num_bands) ; table of fit parameters
t_stat = fltarr(num_bands)
t_signif = fltarr(num_bands)
f_stat = fltarr(num_bands)
f_signif = fltarr(num_bands)

; holdout 1/3 for test purposes
indices1 = where(indgen(n_elements(mask_nochange)) mod 3,complement=indices2)
mask_test = mask_nochange[indices2]
mask_train = mask_nochange[indices1]
n_test = n_elements(mask_test)
n_train = n_elements(mask_train)
n_nochange = n_elements(mask_nochange)

; export regression points to roi
envi_delete_rois, envi_get_roi_ids()
roi_id = envi_create_roi(color=4, name='RadCal regression pixels', ns=num_cols, nl=num_rows)
xpts = mask_train mod num_cols
ypts = mask_train/num_cols
envi_define_roi, roi_id, /point, xpts=xpts, ypts=ypts

print, 'orthogonal regression using ',n_train,' no-change pixels ...'

for i=0,num_bands-1 do begin
   image1 = envi_get_data(fid=fid1,dims=dims1,pos=pos1[i])
   Y = (image1)[mask_train]
   image2 = envi_get_data(fid=fid2,dims=dims2,pos=pos2[i])
   X = (image2)[mask_train]
   window,10+i,xsize=500,ysize=400,xpos=100+10*i,ypos=100+10*i,title="Regression on Band "+strtrim(pos1[i]+1,2)
   wset,10+i
   plot, X, Y, pSym=3,title='Band '+strtrim(pos1[i]+1,2),xtitle="Target",ytitle="Reference", $
       color=0,background='FFFFFF'XL,xrange=[0,max(X)],yrange=[0,max(Y)]
; orthogonal regression on ax+c
   ortho_regress, transpose(X), transpose(Y), ai, xmi, ymi, sigma_aa, sigma_bb, sigma, rank=rank
   aa[i]=ai
   xm[i]=xmi
   ym[i]=ymi
   x0 = min(X) < 0.0
   x1 = max(X)
   Oplot,[x0,x1],[ym[i]-aa[i]*xm[i]+x0*aa[i],ym[i]-aa[i]*xm[i]+x1*aa[i]],color=0
   table[0,i] = i+1              ; band
   table[1,i] = ym[i]-aa[i]*xm[i]; intercept
   table[2,i] = sigma_aa         ; standard error in intercept
   table[3,i] = aa[i]            ; slope
   table[4,i] = sigma_bb         ; standard error in slope
   table[5,i] = correlate(X,Y)   ; correlation
   table[6,i] = sigma            ; RMSE
   table[7,i] = rank[0]          ; Spearman rank-order correlation rho
   table[8,i] = rank[1]          ; significance of difference from zero
; renormalize
   image3[i,*]=ym[i]-aa[i]*xm[i]+image2*aa[i]
endfor
print, 'band, a, sigma_a, b, sigma_b, r, rmse, rho,  s'
print, table

; determine means and standard deviations
; using holdout test pixels only
reference  = fltarr(num_bands,n_test)
target     = fltarr(num_bands,n_test)
normalized = fltarr(num_bands,n_test)
for i=0,num_bands-1 do begin
    reference[i,*] = (envi_get_data(fid=fid1,dims=dims1,pos=pos1[i]))[mask_test]
    target[i,*]    = (envi_get_data(fid=fid2,dims=dims2,pos=pos2[i]))[mask_test]
    normalized[i,*]= (image3[i,*])[mask_test]
    mrs[i]=mean(reference[i,*])
    mts[i]=mean(target[i,*])
    mns[i]=mean(normalized[i,*])
    vrs[i]=variance(reference[i,*])
    vts[i]=variance(target[i,*])
    vns[i]=variance(normalized[i,*])
    tm = tm_test(reference[i,*],normalized[i,*],/paired)
    t_stat[i] = tm[0]
    t_signif[i] = tm[1]
    fv = fv_test(reference[i,*],normalized[i,*])
    f_stat[i] = fv[0]
    f_signif[i] = fv[1]
endfor
; determine covariance matrices and do Wishart test
Sr = correlate(reference,/double,/covariance)
Sn = correlate(normalized,/double,/covariance)
print,"Comparison Statistics using",n_test," test pixels"
print,"Means"
print,"target         ",mts
print,"reference      ",mrs
print,"normalized     ",mns
print,"t-statistic    ",t_stat
print,"p-value        ",t_signif
print,"Variances"
print,"target         ",vts
print,"reference      ",vrs
print,"normalized     ",vns
print,"F-statistic    ",f_stat
print,"p-value        ",f_signif

envi_info_wid,[['RADCAL statistics: '+strtrim(n_train,2)+' training and '+strtrim(n_test,2)+' test pixels'], $
               ['   band      intercept  sigma      slope    sigma        r       RMSE        rho      s'], $
                    string( table, format='(9F10.4)'), $
            'Means',   ['target     '+string(mts,format='('+strtrim(num_bands)+'F10.4)')], $
                       ['reference  '+string(mrs,format='('+strtrim(num_bands)+'F10.4)')], $
                       ['normalized '+string(mns,format='('+strtrim(num_bands)+'F10.4)')], $
                       ['t-stat     '+string(t_stat,format='('+strtrim(num_bands)+'F10.4)')], $
                       ['p-value    '+string(t_signif,format='('+strtrim(num_bands)+'F10.4)')], $
           'Variances',['target     '+string(vts,format='('+strtrim(num_bands)+'F14.6)')], $
                       ['reference  '+string(vrs,format='('+strtrim(num_bands)+'F14.6)')], $
                       ['normalized '+string(vns,format='('+strtrim(num_bands)+'F14.6)')], $
                       ['F-stat     '+string(f_stat,format='('+strtrim(num_bands)+'F14.6)')], $
                       ['p-value    '+string(f_signif,format='('+strtrim(num_bands)+'F14.6)')]], $
          title='RADCAL statistics'

; save to memory or disk
out_array = fltarr(num_cols,num_rows,num_bands)
for i = 0,num_bands-1 do $
    out_array[*,*,i] = reform(image3[i,*],num_cols,num_rows,/overwrite)
if (result.accept eq 0) then begin
  print, 'output cancelled'
end else if (result.outf.in_memory eq 1) then begin
   envi_enter_data, out_array, bnames=bnames2,map_info=map_info2,wl=wl2, xstart=xstart2+dims2[1],ystart=ystart2+dims2[3]
   print, 'result written to memory'
end else begin
   openw, unit, result.outf.name, /get_lun
   for i=0,num_bands-1 do writeu, unit, out_array[*,*,i]
   envi_setup_head ,fname=result.outf.name, ns=num_cols, $
        map_info=map_info2, $
        xstart=xstart2+dims2[1], ystart=ystart2+dims2[3], $
        nl=num_rows, nb=num_bands, $
        data_type=4, interleave=0, /write, $
        bnames=bnames2, $
        wl=wl2, $
        descrip='RadCal normalized'
   print, 'file created ', result.outf.name
   close, unit
endelse

; process another file (e.g. a full scene), writing to disk only

answer=dialog_message('Do you want to use present results to calibrate another file?', /question)
if answer eq 'Yes' then begin
   envi_select, title='Choose file to be calibrated', fid=fid, dims=dims,pos=pos
   if (fid eq -1) then begin
      print,'cancelled'
      return
   endif
   if n_elements(pos) ne num_bands then begin
      print, 'incorrect spectral dimension'
      Message, 'Spectral dimensions are different. Aborting.'
      return
   endif
   envi_file_query, fid, fname=fname, bnames=bnames, interleave=interleave, $
                    xstart=xstart, ystart=ystart, wl=wl
   bnames = bnames[pos]
   wl = wl[pos]                    
   map_info = envi_get_map_info(fid=fid)
   bnames = 'renorm('+bnames+')'
   print,'calibrating file '+fname
   if n_elements(pos1) ne n_elements(pos2) then begin
      print, 'Spectral dimension of file does not match. Aborting.'
      Message, 'Dimensions are different. Aborting.'
      return
   endif
   num_cols = dims[2]-dims[1]+1
   num_rows = dims[4]-dims[3]+1
   base = widget_auto_base(title='Choose output file name')
   sb = widget_base(base, /row, /frame)
   wp = widget_outf(sb, uvalue='outf', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
      print, 'output cancelled'
      return
   endif
   outfilename = result.outf
   openw, unit, outfilename, /get_lun
   progressbar = Obj_New('cgprogressbar',$
                 title='IR-MAD',xsize=300,ysize=20, /cancel)                 
   progressbar->start
; get tile id
   tile_id = envi_init_tile(fid,pos,num_tiles=num_tiles, $
       interleave=interleave,xs=dims[1],xe=dims[2],ys=dims[3],ye=dims[4])
; start tiling
   ym1 = fltarr(num_bands,num_cols)
   xm1 = fltarr(num_bands,num_cols)
   aa1 = fltarr(num_bands,num_cols)
   for i=0,num_cols-1 do begin
      ym1[*,i] = ym
      xm1[*,i] = xm
      aa1[*,i] = aa
   endfor
   for tile_index=0L,num_tiles-1 do begin
      tile = envi_get_tile(tile_id,tile_index,band_index=k)
      case interleave of
        0: image1 = ym[k] - aa[k]*xm[k] + aa[k]*tile
        1: image1 = transpose(ym1-aa1*xm1) + transpose(aa1)*tile
        2: image1 = ym1 - aa1*xm1 + aa1*tile
      endcase
      if progressbar->CheckCancel() then begin
         print,'Calibration aborted'
         free_lun, unit
         progressbar->Destroy
         return
      endif
      progressbar->Update,tile_index*100/num_tiles
      writeu,unit,image1
   endfor
; tidy up
   free_lun, unit
   envi_tile_done, tile_id
   progressbar->Destroy
   envi_setup_head,fname=outfilename, ns=num_cols, $
                   nl=num_rows, nb=num_bands, $
                   data_type=4, $
                   file_type=0,map_info=map_info, $
                   interleave=interleave, /write, $
                   xstart=xstart+dims[1], $
                   ystart=ystart+dims[3], $
                   wl=wl,$
                   bnames=bnames
   print, 'calibrated file created ', outfilename
endif

end