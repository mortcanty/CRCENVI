 ;docformat = 'rst'
; c_correction_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.


PRO c_correction_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'C-Correction', $
      REF_VALUE = '3D SurfaceView', $
      EVENT_PRO = 'c_correction_run', $
      UVALUE = 'C_CORRECTION',$
      POSITION = 'after', $
      /SEPARATOR
END

;+
; :Description:
;       c-correction for solar
;       illumination in rough terrain.
;        dem projection (and rotation)
;        must be same as that of image
;        to be corrected.
;       
;       Ref: P. M. Teillet et al. (1982), Canadian Journal of
;       Remote Sensing 8(2) 84-106    
; :Params:
;      event:  in, required 
;         if called from the ENVI menu                  
; :Uses:
;      ENVI
; :Author:
;      Mort Canty (2013)
;      m.canty@fz-juelich.de        
;-
pro c_correction_run, event

COMPILE_OPT IDL2

print, '---------------------------------'
print, 'Solar illumination c-correction'
print, systime(0)
print, '---------------------------------'

catch, theError
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
   return
endif

envi_select, title='Choose image to correct', fid=fid_image, pos=pos_image, dims=dims_image, $
                                              /mask,m_fid=fid_mask,m_pos=pos_mask
if (fid_image eq -1) then begin
   print, 'cancelled'
   return
endif

envi_file_query, fid_image, fname=fname, bnames=bnames, xstart=xstart, ystart=ystart
num_bands =  n_elements(pos_image)
num_cols = dims_image[2]-dims_image[1]+1
num_rows = dims_image[4]-dims_image[3]+1
num_pixels = num_cols*num_rows

; tie point
map_info = envi_get_map_info(fid=fid_image)
envi_convert_file_coordinates, fid_image, dims_image[1], dims_image[3], e, n, /to_map
map_info.mc= [0D,0D,e,n]

if (fid_mask eq -1) then  begin
   mask = lindgen(num_pixels)
   nomask = -1
   ncomplement = 0
end else mask=where(envi_get_data(fid=fid_mask, dims=dims_image, pos=pos_mask),complement=nomask,ncomplement=ncomplement)

envi_select, title='Choose digital elevation file', fid=fid_dem, dims=dims, pos=pos, /band
if (fid_dem eq -1) then begin
   print, 'cancelled'
   return
endif

map_info_dem = envi_get_map_info(fid=fid_dem)
; check for same projection
if (map_info.proj.name ne map_info_dem.proj.name) or (round(map_info.rotation) ne round(map_info_dem.rotation)) then $
   message, 'map projections do not match' 

; get resize factor for dem
rfact = map_info.ps/map_info_dem.ps

base = widget_auto_base(title='C-correction parameters')
list = ['solar_elevation ', 'solar_azimuth', 'topo kernel']
vals = [0.0,0.0,9.0]
we = widget_edit(base,  list=list, uvalue='edit', $
      vals=vals, field= 2, dt=4, /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print,'cancelled'
   return
endif
print, 'solar elevation(deg) ', result.edit[0]
print, 'solar azimuth(deg) ', result.edit[1]
print, 'topographic kernel ', result.edit[2]
theta_z = (90-result.edit[0])*!DTOR ; solar zenith angle
phi_a = result.edit[1]*!DTOR        ; solar azimuth angle
kernel = fix(result.edit[2])        ; kernel

; output destination
base = widget_auto_base(title='Output corrected image')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
    print, 'cancelled'
    return
endif

; find upper left position of image within DEM
envi_convert_file_coordinates, fid_dem, X_ul, Y_ul, e, n

; cutout the corresponding spatial subset of the DEM and place in memory
dims1=[-1L, fix(X_ul), fix(X_ul)+num_cols*rfact[0]-1, fix(Y_ul), fix(Y_ul)+num_rows*rfact[1]-1]
dem = envi_get_data(fid=fid_dem, dims=dims1, pos=0)
envi_enter_data, dem, r_fid=r_fid

; resize the DEM and place in memory
dims1 = [-1L, 0, num_cols*rfact[0]-1, 0, num_rows*rfact[1]-1]
envi_doit, 'resize_doit', fid=r_fid, dims=dims1, pos=0, /in_memory, interp=0, r_fid=r_fid1, rfact=rfact
envi_file_mng, id=r_fid, /remove

;determine slope and aspect from DEM
dims1 = [-1L, 0, num_cols-1, 0, num_rows-1]
envi_doit, 'topo_doit', bptr=[0,1], dims=dims1, pos=0, kernel=kernel, $
          fid=r_fid1, /in_memory, pixel_size=map_info.ps, r_fid=r_fid2
envi_file_mng, id=r_fid1, /remove

; calculate local solar elevation
theta_p = envi_get_data(fid=r_fid2, dims=dims1, pos=[0])*!DTOR
phi_o   = envi_get_data(fid=r_fid2, dims=dims1, pos=[1])*!DTOR
cos_gamma_i = cos(theta_p)*cos(theta_z)+sin(theta_p)*sin(theta_z)*cos(phi_a-phi_o)
envi_file_mng, id=r_fid2, /remove

; save cos_gamma_i image to ENVI
envi_enter_data, cos_gamma_i+0.0, map_info=map_info, descrip='cos(gamma_i)'

; apply the c-correction
image=fltarr(num_pixels, num_bands)
widget_control, /hourglass
for k=0, num_bands-1 do begin
   band = envi_get_data(fid=fid_image, dims=dims_image, pos=pos_image[k])
   m=( regress(cos_gamma_i[mask], band[mask], const=b, correlation=r, /double) )[0]
   print,'band '+strtrim(k+1,2)+': m = '+strtrim(m,2)+'  b = '+strtrim(b,2)+'  r = '+strtrim(r[0],2)  
   image[mask,k] = band[mask]*(cos(theta_z) + b/m)/(cos_gamma_i[mask] + b/m)
   if ncomplement gt 0 then image[nomask,k] = band[nomask]
endfor

image = reform(image,num_cols,num_rows,num_bands)

if (result.outf.in_memory eq 1) then begin
   envi_enter_data, image, $
        bnames=bnames[pos], $
        map_info=map_info, $
        xstart=xstart+dims_image[1], ystart=ystart+dims_image[3], $
        descrip='C-Correction to image:'+file_basename(result.outf.name)
   print, 'Result written to memory'
endif else begin
   openw, unit, result.outf.name, /get_lun
   for i=0,num_bands-1 do writeu, unit, image[*,*,i]
   envi_setup_head ,fname=result.outf.name, ns=num_cols, $
        nl=num_rows, nb=num_bands, $
        data_type=4, interleave=0, /write, $
        bnames=bnames[pos], $
        map_info=map_info, $
        xstart=xstart+dims_image[1], ystart=ystart+dims_image[3], $
        descrip='C-Correction to image:'+file_basename(result.outf.name)
   print, 'File created ', result.outf.name
   free_lun, unit
endelse

end