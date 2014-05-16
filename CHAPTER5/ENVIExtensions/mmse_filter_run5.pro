; docformat = 'rst'
; mmse_filter_run5.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
    
pro mmse_filter_run5_extensions_init
   compile_opt idl2
   e = envi(/current)
   e.addextension, 'MMSE polSAR filter', 'mmse_filter_run5', path='Chapter5'
end

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
;    Lee MMSE adaptive filtering 
;    for polSAR covariance images
;    Lee et al. (1999) IEEE TGARS 37(5), 2363-2373
;    Oliver and Quegan (2004) Understanding SAR Images, Scitech 
;    Should be valid for all matrix elements
; :Params:
;    NONE  
; :KEYWORDS:
;    NONE    
; :Uses:
;    COYOTE, ENVI    
; :Author:
;       Mort Canty (2013)        
;-
pro mmse_filter_run5

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
    
   print, '--------------------------------------'
   print, 'MMSE filter of polSAR covariance image'
   print, systime(0)
   print, '--------------------------------------'
    

   
   templates = intarr(7,7,8)
   for j = 0,6 do begin
       templates[0:2,j,0] = 1
       templates[0,1,1] = 1
       templates[0:1,2,1] = 1
       templates[0:2,3,1] = 1
       templates[0:3,4,1] = 1
       templates[0:4,5,1] = 1
       templates[0:5,6,1] = 1
   endfor    
   templates[*,*,6] = transpose(templates[*,*,0])
   templates[*,*,2] = reverse(templates[*,*,6],2)
   templates[*,*,3] = reverse(templates[*,*,1],1)
   templates[*,*,4] = reverse(templates[*,*,0],1)
   templates[*,*,5] = reverse(templates[*,*,3],2)
   templates[*,*,7] = reverse(templates[*,*,5],1)   
   templates = reform(templates,49,8)
   tmp = intarr(21,8)
   for i = 0,7 do tmp[*,i] = where(templates[*,i]) 
   templates = tmp
   
   edges = intarr(3,3,4)
   edges[*,*,0] = [[-1,0,1],[-1,0,1],[-1,0,1]]
   edges[*,*,1] = [[0,1,1],[-1,0,1],[-1,-1,0]]
   edges[*,*,2] = [[1,1,1],[0,0,0],[-1,-1,-1]]
   edges[*,*,3] = [[1,1,0],[1,0,-1],[0,-1,-1]]   
   
   es = fltarr(4)
   
   ; input multi-look averaged covariance image
   envi_select, title='Choose (spatial subset of) covariance matrix image', $
     fid=fid, dims=dims, pos=pos
   if (fid eq -1) then begin
     print, 'cancelled'
     return
   end
   envi_file_query, fid, fname=fname, bnames=bnames,xstart=xstart, ystart=ystart
   cols = dims[2]-dims[1]+1
   rows = dims[4]-dims[3]+1
   bands = n_elements(pos)
   ; bands = 9: quad pol
   ; bands = 4  dual pol
   ; bands = 1  single pol
   
   ; tie point
   map_info = envi_get_map_info(fid=fid)
   envi_convert_file_coordinates, fid, dims[1], dims[3], e, n, /to_map
   map_info.mc= [0D,0D,e,n]
   
   ; number of looks
   base = widget_auto_base(title='Number of looks')
   we = widget_param(base, dt=4, field=5, floor=1, $
     default=5, uvalue='param', xsize=32, /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then m=5 else m=result.param
   
   ; output destination
   base = widget_auto_base(title='Output filtered image')
   sb = widget_base(base, /row, /frame)
   wp = widget_outfm(sb, uvalue='outf', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
     print, 'cancelled'
     return
   endif
   
   ; construct span image
   span = envi_get_data(fid=fid,dims=dims,pos=pos[0])
   if bands eq 9 then begin
     span = span + envi_get_data(fid=fid,dims=dims,pos=pos[5])
     span = span + envi_get_data(fid=fid,dims=dims,pos=pos[8])
   end else if bands eq 4 then span = span + envi_get_data(fid=fid,dims=dims,pos=pos[3])
   
   ; get the filter weights from span image
   print, 'Determining filter weight image ...'
   b = fltarr(cols,rows) + 1.0
   edge_idx = intarr(cols,rows)
   start = systime(2)
   progressbar = Obj_New('cgprogressbar', /cancel, title='MMSE ...')
   progressbar->start
   for j = 3L,rows-4 do begin
     if progressbar->CheckCancel() then begin
       print,'ENL interrupted'
       progressbar->Destroy
       return
     endif
     get_windex, j, cols, windex
     for i = 3L,cols-4 do begin
       wind = reform(span[windex],7,7)
       ;       3x3 compression
       w = congrid(wind,3,3,/center,cubic = -0.5)
       ;       get appropriate edge mask
       for p=0,3 do es[p] = total(edges[*,*,p]*w)
       void = max(es,idx)
       case idx of
         0: if abs(w[1,1]-w[0,1]) lt abs(w[1,1]-w[2,1]) then $
           edge_idx[i,j] = 0 else $
           edge_idx[i,j] = 4
         1: if abs(w[1,1]-w[0,2]) lt abs(w[1,1]-w[2,0]) then $
           edge_idx[i,j] = 1 else $
           edge_idx[i,j] = 5
         2: if abs(w[1,1]-w[1,0]) lt abs(w[1,1]-w[1,2]) then $
           edge_idx[i,j] = 6 else $
           edge_idx[i,j] = 2
         3: if abs(w[1,1]-w[0,0]) lt abs(w[1,1]-w[2,2]) then $
           edge_idx[i,j] = 7 else $
           edge_idx[i,j] = 3
       endcase
       edge = templates[*,edge_idx[i,j]]
       wind = wind[edge]
       gbar = mean(wind)
       varg = variance(wind)
       b[i,j] = (1.0 - gbar^2/(varg*m))/(1.0+1.0/m) > 0.0
       windex++
     endfor
     progressbar->Update,j*100/rows
   endfor
   progressbar->destroy
   ; filter the image
   outimage = fltarr( cols,rows,bands )
   envi_file_query,fid,dims=dims
   print, 'Filtering covariance matrix elements ...'
   for k = 0,bands-1 do begin
     print, 'band '+strtrim(k+1,2)
     inband = envi_get_data(/complex,fid=fid,dims=dims,pos=k)
     gbar = inband*0.0
     ;    get window means
     for j = 3L,rows-4 do begin
       get_windex, j, cols, windex
       for i = 3L,cols-4 do begin
         wind = inband[windex]
         edge = templates[*,edge_idx[i,j]]
         wind = wind[edge]
         gbar[i,j] = mean(wind)
         windex++
       endfor
     endfor
     ;    apply adaptive filter
     outimage[*,*,k] = gbar + b*(inband-gbar)
   endfor
   
   if (result.outf.in_memory eq 1) then begin
     envi_enter_data, outimage, $
       bnames=bnames[pos], $
       map_info=map_info, $
       xstart=xstart+dims[1], ystart=ystart+dims[3], $
       descrip='MMSE filter:'+file_basename(result.outf.name)
     print, 'Result written to memory'
   endif else begin
     openw, unit, result.outf.name, /get_lun
     for i=0,bands-1 do writeu, unit, outimage[*,*,i]
     envi_setup_head ,fname=result.outf.name, ns=cols, $
       nl=rows, nb=bands, $
       data_type=4, interleave=0, /write, $
       bnames=bnames[pos], $
       map_info=map_info, $
       xstart=xstart+dims[1], ystart=ystart+dims[3], $
       descrip='MMSE filter:'+file_basename(result.outf.name)
     print, 'File created ', result.outf.name
     free_lun, unit
   endelse
   envi_file_mng, /remove, id=fid
   print, 'done, elapsed time: '+strtrim(systime(2)-start,2)+' sec'

end