; docformat = 'rst'
; bslfcpnorm_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO bslfcpnorm_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Scatterplot Normalization', $
      REF_VALUE = 'Change Detection Statistics', $
      EVENT_PRO = 'bslfcpnorm_run', $
      UVALUE = 'BSLFCPNORM',$
      POSITION = 'after', $
      /SEPARATOR
END

;+
; :Description:
;       ENVI extension for BSL FCP normalization::
; :Params:
;       event:  in, optional 
;          required if called from ENVI  
; :Uses:
;       ENVI               
; :Author:
;       Mort Canty (2013)
;-
pro bslfcpnorm_run, event

COMPILE_OPT IDL2 

    print, '---------------------------------'
    print, 'Scatterplot normalization'
    print, systime(0)
    print, '---------------------------------'
    
    catch, theError
    if theError ne 0 then begin
       void = Dialog_Message(!Error_State.Msg, /error)
       return
    endif
    
    envi_select, title='Choose RED and NIR bands of reference image', fid=fid1, dims=dims1,pos=pos1
        if (fid1 eq -1) then begin
          print,'cancelled'
          return
        endif
    envi_file_query, fid1, fname=fname1, ns=ns, nl=nl, nb=nb
    print, 'Reference: ',fname1
    envi_select, title='Choose RED and NIR bands of target image', fid=fid2, dims=dims2,pos=pos2
        if (fid2 eq -1) then begin
           print,'cancelled'
           return
        endif
    envi_file_query, fid2, fname=fname2, bnames=bnames2, ns=ns, nl=nl, nb=nb, xstart=xstart2, ystart=ystart2, wl=wl2
    print, 'Target:    ',fname2
    print, '---------------------------------'
    if ((n_elements(pos1) ne 2) or (n_elements(pos2) ne 2)) then $
       Message, 'Number of selected bands must be 2. Aborting.'
    num_cols1 = dims1[2]-dims1[1]+1
    num_rows1 = dims1[4]-dims1[3]+1
    num_pixels1 = (num_cols1*num_rows1)
    bnames2 = 'renorm('+bnames2[pos2]+')'
    wl2 = wl2[pos2]
    num_cols2 = dims2[2]-dims2[1]+1
    num_rows2 = dims2[4]-dims2[3]+1
    num_pixels2 = (num_cols2*num_rows2)
    
; output destination
    base = widget_auto_base(title='Save normalized image')
    sb = widget_base(base, /row, /frame)
    wp = widget_outfm(sb, uvalue='outf', /auto)
    result = auto_wid_mng(base)
    
; tie point
    map_info2 = envi_get_map_info(fid=fid2)
    envi_convert_file_coordinates, fid2, dims2[1], dims2[3], e, n, /to_map
    map_info2.mc= [0D,0D,e,n]
    
; read in images
    r_image = fltarr(num_pixels2,2)
    t_image = fltarr(num_pixels2,2)
    for i=0,1 do begin
       r_image[*,i] = envi_get_data(fid=fid1,dims=dims1,pos=pos1[i])
       t_image[*,i] = envi_get_data(fid=fid2,dims=dims2,pos=pos2[i])
    endfor
    
; BSL and FCP for reference image
    window,11,xsize=500,ysize=400,xpos=100,ypos=100,title="REFERENCE"
    wset,11
    plot, r_image[0:*:10,0],r_image[0:*:10,1] , pSym=3,xtitle="RED",ytitle="NIR", $
           color=0,background='FFFFFF'XL,xrange=[0,250],yrange=[0,250]
    idx=where(r_image[*,0] gt 0 and r_image[*,1] le 255)
    ratio=r_image[idx,1]/r_image[idx,0]
    nbins=1000
    binsize=(max(ratio) - min(ratio))/nbins
    hist=histogram(ratio,binsize=binsize, $
                               reverse_indices=ri)
    i = 1
; 0.01 quantile
    thresh = num_pixels1/100L
    while ri[i]-ri[0] lt thresh do i++
    idx1 = idx[ ri[ri[0]:ri[i-1]-1] ]
    X = r_image[idx1,0]
    Y = r_image[idx1,1]
    ortho_regress,transpose(X),transpose(Y),A_R,xm,ym
    B_R = ym-A_R*xm       
    oplot, [0,255],[B_R,255*A_R+B_R],color='0000FF'XL     
; 0.999 quantile
    thresh = num_pixels1/1000L
    i= 999  
    while ri[1000]-ri[i] lt thresh do i--
    idx2 = idx[ ri[ri[i]:*] ]
    FCP = mean( r_image[idx2,*],dimension=1 )
    X_R = FCP[0]
    Y_R = FCP[1]
    oplot, [X_R],[Y_R],psym=1,color='0000FF'XL
    
; BSL and FCP for target image
    window,12,xsize=500,ysize=400,xpos=200,ypos=200,title="TARGET"
    wset,12
    plot, t_image[0:*:10,0],t_image[0:*:10,1] , pSym=3,xtitle="RED",ytitle="NIR", $
           color=0,background='FFFFFF'XL,xrange=[0,250],yrange=[0,250]
    idx = where(t_image[*,0] gt 0 and t_image[*,1] le 255)
    ratio = t_image[idx,1]/t_image[idx,0]
    nbins = 1000
    binsize = (max(ratio) - min(ratio))/nbins
    hist = histogram(ratio,binsize=binsize,reverse_indices=ri)
    i = 1
; 0.01 quantile
    thresh = num_pixels1/100L
    while ri[i]-ri[0] lt thresh do i++
    idx1 = idx[ ri[ri[0]:ri[i-1]-1] ]
    X = t_image[idx1,0]
    Y = t_image[idx1,1]
    ortho_regress, transpose(X), transpose(Y), A_T, xm, ym
    B_T = ym-A_T*xm       
    oplot, [0,255],[B_T,255*A_T+B_T],color='0000FF'XL     
; 0.999 quantile
    thresh = num_pixels1/1000L
    i= 999  
    while ri[1000]-ri[i] lt thresh do i--
    idx2 = idx[ ri[ri[i]:*] ]
    FCP = mean( t_image[idx2,*],dimension=1 )
    X_T = FCP[0]
    Y_T = FCP[1]
    oplot, [X_T],[Y_T],psym=1,color='0000FF'XL
    
; normalize
    L_R = Y_R - A_R*X_R+B_R
    L_T = Y_T - A_T*X_T+B_T
    Y = Y_R - (Y_T - t_image[*,1])*L_R/L_T
    X = X_R + (t_image[*,0] - X_T)*(L_R*A_T)/(L_T*A_R)
    
    window,13,xsize=500,ysize=400,xpos=300,ypos=300,title="NORMALIZED TARGET"
    wset,13
    plot, X[0:*:10],Y[0:*:10] , pSym=3,xtitle="RED",ytitle="NIR", $
           color=0,background='FFFFFF'XL,xrange=[0,250],yrange=[0,250]
    
; write result
    out_array = fltarr(num_cols2,num_rows2,2)
    out_array[*,*,0] = reform(X,num_cols2,num_rows2,/overwrite)
    out_array[*,*,1] = reform(Y,num_cols2,num_rows2,/overwrite)
    if (result.accept eq 0) then begin
      print, 'output cancelled'
    end else if (result.outf.in_memory eq 1) then begin
       envi_enter_data, out_array, bnames=bnames,map_info=map_info2,wl=wl2, xstart=xstart2+dims2[1],ystart=ystart2+dims2[3]
       print, 'result written to memory'
    end else begin
       openw, unit, result.outf.name, /get_lun
       for i=0,1 do writeu, unit, out_array[*,*,i]
       envi_setup_head ,fname=result.outf.name, ns=num_cols2, $
            map_info=map_info2, $
            xstart=xstart2+dims2[1], ystart=ystart2+dims2[3], $
            nl=num_rows2, nb=2, $
            data_type=4, interleave=0, /write, /open, $
            bnames=bnames2, $
            wl=wl2, $
            descrip='RadCal normalized'
       print, 'file created ', result.outf.name
       close, unit
    endelse

End