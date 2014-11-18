; docformat = 'rst'
; som_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO som_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Self Organizing Map', $
      REF_VALUE = 'K-Means', $
      EVENT_PRO = 'som_run', $
      UVALUE = 'SOM',$
      POSITION = 'after', $
      /SEPARATOR
END

function winner, g, W
   sz = size(W)
   GG = W*0.0
   for i=0,sz[2]-1 do replicate_inplace, GG, g[i], 1, [0,i]
   D = total((W - GG)^2,2)
   return,(where(min(D) eq D))[0]
end

function dsquare, m, ks, k
  return, float( (((ks-1) mod m) - ((k-1) mod m))^2 + (((ks-1)/m - (k-1)/m) mod m)^2 + ((ks-1)/m^2 - (k-1)/m^2)^2 )
end

;+
; :Description:
;       ENVI extension for Kohonen Self 
;       Organizing Map::
;           Kohonen, T. (1989). Self-Organization 
;           and Associative Memory. Springer.
; :Params:
;       event:  in, optional 
;          required if called from ENVI                 
; :Uses:
;       ENVI::
;       COYOTE            
; :Author:
;       Mort Canty (2013)      
pro SOM_run, event

COMPILE_OPT STRICTARR

print, '---------------------------------'
print, 'SOM Clustering'
print, systime(0)
print, '---------------------------------'

catch, theError
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
   return
endif

envi_select, title='Choose multispectral image', $
             fid=fid, dims=dims,pos=pos
if (fid eq -1) then begin
   print, 'cancelled'
   return
endif

envi_file_query, fid, fname=fname, xstart=xstart, ystart=ystart
print, 'Selected image: ',fname

; tie point
map_info = envi_get_map_info(fid=fid)
envi_convert_file_coordinates, fid, dims[1], dims[3], e, n, /to_map
map_info.mc[2:3]= [e,n]

base = widget_auto_base(title='Cube side dimension')
wg = widget_sslider(base, title='Dimension', min=4, max=10, $
  value=6, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then return

base = widget_auto_base(title='SOM Output')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result1 = auto_wid_mng(base)
if (result1.accept eq 0) then begin
  print, 'Output cancelled'
  return
endif

K = long(result.slide)
print, 'Cube side dimension', K

num_cols = dims[2]-dims[1]+1L
num_rows = dims[4]-dims[3]+1L
num_pixels = (num_cols*num_rows)
num_bands = n_elements(pos)
if num_bands lt 3 then message, 'Dimnesion of image must be at least 3'

widget_control, /hourglass

; random sample of 10000 pixels
indices = randomu(seed, 10000, /long) mod num_pixels
samples = fltarr(n_elements(indices), num_bands)
image = fltarr(num_pixels, num_bands)
for i=0,num_bands-1 do begin
    temp= envi_get_data(fid=fid,dims=dims,pos=pos[i])
    image[*,i] = temp
    samples[*,i]=temp[indices]
endfor

; initialize synaptic weights with training data
indices = randomu(seed,K^3,/long) mod K^3
W = samples[indices,*]

; train the network
print, 'training...'
etamax = 1.0
etamin = 0.001
sigmamax = K/2.0
sigmamin = 0.5
progressbar = Obj_New('cgprogressbar',$
              title='SOM training...',xsize=250,ysize=20, /cancel)
progressbar->start
for i=0L,9999 do begin
   if progressbar->CheckCancel() then begin
      print,'training aborted'
      goto, done
   endif
   progressbar->Update,(i*100)/10000
   g = samples[i,*]
   kstar = winner(g, W)
   eta = etamax*(etamin/etamax)^(i/10000.0)
   sigma = sigmamax*(sigmamin/sigmamax)^(i/10000.0)
   for j=0L,K^3-1 do begin
      d2 = dsquare(K,kstar,j)
      lambda = exp(-d2/(2*sigma^2))
      W[j,*] = W[j,*] + eta*lambda*(g-W[j,*])
   endfor
endfor

done:
progressbar->destroy

; display the neurons in 3-D feature space
symbols = objarr(K^3)
for i=0L,K^3-1 do begin
   red = i mod K
   green = i/K mod K
   blue = i/K^2
   color = byte([fix(red*255/(K-1)),fix(green*255/(K-1)),fix(blue*255/(K-1))])
   oOrb = OBJ_NEW('orb', color=color)
   oOrb->Scale, 1.0, 1.0, 1.0
   symbols[i] = OBJ_NEW('IDLgrSymbol', oOrb)
end
xplot3d, W[*,0],W[*,1],W[*,2],linestyle=6,symbol=symbols

; cluster the image
progressbar = Obj_New('cgprogressbar', $
              title= 'SOM classifying...',xsize=250,ysize=20)
progressbar->start
class_image = bytarr(num_pixels,3)
for i=0L,num_pixels-1 do begin
   if i mod 1000 eq 0 then begin
      if progressbar->CheckCancel() then begin
         print,'classifying aborted'
         goto, done1
      endif
      incr = (i*100)/num_pixels
      progressbar->Update,incr
   endif
   x = image[i,*]
   kstar = winner(x,W)
   red = kstar mod K
   class_image[i,0] = fix(red*255/(K-1))
   green = kstar/K mod K
   class_image[i,1] = fix(green*255/(K-1))
   blue = kstar/K^2
   class_image[i,2] = fix(blue*255/(K-1))
endfor

done1:
progressbar->destroy

; output to memory or file
if (result1.outf.in_memory eq 1) then begin
   envi_enter_data, reform(class_image,num_cols,num_rows,3), $
                    file_type=1, $
                    map_info=map_info, $
                    bnames=['SOM('+fname+')'],$
                    xstart=xstart+dims[1], $
                    ystart=ystart+dims[3]
   print, 'Result written to memory'
endif else begin
   openw, unit, result1.outf.name, /get_lun
   writeu, unit, reform(class_image,num_cols,num_rows,3)
   envi_setup_head ,fname=result1.outf.name, ns=num_cols, nl=num_rows, nb=3, $
       data_type=1, $
       interleave=0, $
       file_type=0, $
       map_info=map_info, $
       xstart=xstart+dims[1], $
       ystart=ystart+dims[3], $
       descrip='SOM clustering of ' + result1.outf.name, $
       /write
   print, 'File created ', result1.outf.name
   free_lun, unit
endelse

end