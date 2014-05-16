; docformat = 'rst'
; gamma_filter_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
    
PRO gamma_filter_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Gamma MAP polSAR filter', $
      REF_VALUE = 'AIRSAR Scattering Classification', $
      EVENT_PRO = 'gamma_filter_run', $
      UVALUE = 'GAMMA',$
      POSITION = 'after'
END

pro get_windex, j, cols, windex
;  7x7 window, centered at (3,j)
    windex = lonarr(49)
    windex[0:6]   = (j-3)*cols + [0,1,2,3,4,5,6]
    windex[7:13]  = (j-2)*cols + [0,1,2,3,4,5,6]
    windex[14:20] = (j-1)*cols + [0,1,2,3,4,5,6]
    windex[21:27] = (j)*cols   + [0,1,2,3,4,5,6] 
    windex[28:34] = (j+1)*cols + [0,1,2,3,4,5,6] 
    windex[35:41] = (j+2)*cols + [0,1,2,3,4,5,6]
    windex[42:48] = (j+3)*cols + [0,1,2,3,4,5,6]   
end

function gamma_filter, k, result=result
   compile_opt strictarr
   common gfilt, templates, edges, inimage, outimage, rows, cols, outbnames, m
   result = inimage[*,*,k]
   arr = outimage[*,*,k]
   es = fltarr(4)
   print, 'filtering band '+outbnames[k]
   progressbar = Obj_New('cgprogressbar', /cancel, title='Gamma MAP...')  
   progressbar->start      
   for j = 3L,rows-4 do begin
      if progressbar->CheckCancel() then begin
         print,'Filter interrupted'
         progressbar->Destroy
         return, 0
      endif   
      get_windex, j, cols, windex
      for i = 3L,cols-4 do begin
;       central pixel         
         g = inimage[i,j,k]      
;       content of 7x7 window         
         wind = reform(arr[windex],7,7)  
;       3x3 compression         
         w = congrid(wind,3,3,/center,cubic = -0.5)
;       get appropriate edge mask
         for p=0,3 do es[p] = total(edges[*,*,p]*w)
         void = max(es,idx)
         case idx of
            0: if abs(w[1,1]-w[0,1]) lt abs(w[1,1]-w[2,1]) then $
                    edge = templates[*,0] else $
                    edge = templates[*,4]
            1: if abs(w[1,1]-w[0,2]) lt abs(w[1,1]-w[2,0]) then $
                    edge = templates[*,1] else $
                    edge = templates[*,5]   
            2: if abs(w[1,1]-w[1,0]) lt abs(w[1,1]-w[1,2]) then $
                    edge = templates[*,6] else $
                    edge = templates[*,2]
            3: if abs(w[1,1]-w[0,0]) lt abs(w[1,1]-w[2,2]) then $
                    edge = templates[*,7] else $
                    edge = templates[*,3]
         endcase
         wind = wind[edge]
         var = variance(wind)
         mu = mean(wind)         
         alpha = abs( (1 +1/m)/(var/mu^2 - 1/m) )
         a = mu*(alpha-m-1)
         x = (a+sqrt(4*g*m*alpha*mu+a^2))/(2*alpha)        
         result[i,j] = x
         windex++  
      endfor   
      progressbar->Update,j*100/rows
   endfor  
   progressbar->destroy    
   return, 1
end    
    
;+
; :Description:   
;    (iterated) gamma MAP adaptive filtering for polarized SAR intensity images
;    Ref: Oliver and Quegan (2004) Understanding SAR Images, Scitech 
;    Valid only for diagonal covariance matrix elements
; :Params:
;      event:  in, required 
;         if called from the ENVI Classic menu   
; :KEYWORDS:
;    NONE   
; :Uses:
;    POLSAR__DEFINE, ENVI    
; :Author:
;       Mort Canty (2013)        
;-
pro gamma_filter_run, event 
   compile_opt idl2
   common gfilt, templates, edges, inimage, outimage, rows, cols, outbnames, m
    
   print, '-------------------------------------------'
   print, 'Gamma MAP filter of polSAR covariance image'
   print, systime(0)
   print, '-------------------------------------------'
    
   catch, theError
   if theError ne 0 then begin
      void = Dialog_Message(!Error_State.Msg, /error)
      return
   endif
 
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
   
   print, '-----------------------'
   print, 'Gamma MAP polSAR filter'
   print, systime(0)
   print, '-----------------------'   
    
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
   
; number of iterations   
   base = widget_auto_base(title='Number of iterations')
   we = widget_param(base, dt=1, field=5, floor=1, $
      default=1, uvalue='param', xsize=32, /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then niter=1 else niter=round(result.param)   
   print, 'Number of iterations '+strtrim(niter,2)  
   
   ; output destination
   base = widget_auto_base(title='Output filtered image')
   sb = widget_base(base, /row, /frame)
   wp = widget_outfm(sb, uvalue='outf', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
     print, 'cancelled'
     return
   endif   
   
; process the diagonal matrix elements only 
   if bands eq 9 then outbands = 3 $
   else if bands eq 4 then outbands = 2 $
   else if bands eq 1 then outbands = 1 $
   else begin
     print, 'Incorrect number of bands'
     message, 'Incorrect number of bands'
   endelse
   
   case outbands of
    3: outbnames = ['C11','C22','C33']
    2: outbnames = ['C11','C22']
    1: outbnames = ['C11']
   endcase

   inimage = fltarr(cols,rows,outbands) 
   inimage[*,*,0] = envi_get_data(fid=fid,dims=dims,pos=pos[0])
   if outbands eq 2 then inimage[*,*,1] = envi_get_data(fid=fid,dims=dims,pos=pos[3]) $
   else if outbands eq 3 then begin
      inimage[*,*,1] = envi_get_data(fid=fid,dims=dims,pos=pos[5])
      inimage[*,*,2] = envi_get_data(fid=fid,dims=dims,pos=pos[8])
   endif 
   outimage = inimage + 0.0
   
   itr = 0
   start = systime(2)

   while itr lt niter do begin
      print, 'iteration '+strtrim(itr+1,2)
      for i = 0,outbands-1 do begin 
        if not gamma_filter(i,result=res) then begin
          print, 'interrupted, aborting ...'
          return
        endif
        outimage[*,*,i] = res
      endfor             
      itr++
   endwhile  
   
   if (result.outf.in_memory eq 1) then begin
     envi_enter_data, outimage, $
       bnames=outbnames, $
       map_info=map_info, $
       xstart=xstart+dims[1], ystart=ystart+dims[3], $
       descrip='Gamma MAP filter:'+file_basename(result.outf.name)
     print, 'Result written to memory'
   endif else begin
     openw, unit, result.outf.name, /get_lun
     for i=0,outbands-1 do writeu, unit, outimage[*,*,i]
     envi_setup_head ,fname=result.outf.name, ns=cols, $
       nl=rows, nb=outbands, $
       data_type=4, interleave=0, /write, $
       bnames=outbnames, $
       map_info=map_info, $
       xstart=xstart+dims[1], ystart=ystart+dims[3], $
       descrip='Gamma MAP filter:'+file_basename(result.outf.name)
     print, 'File created ', result.outf.name
     free_lun, unit
   endelse
   envi_file_mng, /remove, id=fid
   print, 'done, elapsed time: '+strtrim(systime(2)-start,2)+' sec'   

end   
   
   