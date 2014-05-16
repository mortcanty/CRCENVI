; docformat = 'rst'
; mad_run_stats.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

;+
; :Description:
;       ENVI extension for Iteratively Re-weighted
;       Multivariate Alteration Detection (IR-MAD)
;       using existing satistics file::
; :Params:
;       event:  in, optional 
;          required if called from ENVI  
; :Uses:
;       ENVI::            
;       COYOTE                 
; :Author:
;       Mort Canty (2009)
;-
pro mad_run_stats, event

COMPILE_OPT STRICTARR 

print, '---------------------------------'
print, 'Multivariate Alteration Detection'
print, systime(0)
print, '---------------------------------'

catch, theError
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
   return
endif

; stats input
fileName = Dialog_Pickfile(Filter=['*.mad'],/Read,Title='Read existing MAD stats file')
if fileName eq '' then begin
    Print, 'Cancelled'
    return
endif     
; read MAD stats file
print, 'reading MAD stats'  
openR, lun, fileName, /get_lun
; read and print the two-line header
aLine = ''
for i=0,1 do begin
   readF, lun, aLine
   print, aLine
endfor
; get number of bands
nbands = 0L 
readF, lun, nbands
means1 = (means2 = (rho = (sigma = fltarr(nbands))))
A = (B = fltarr(nbands,nbands))
readF, lun, means1, means2, A, B, rho, sigma
free_lun, lun

envi_select, title='Choose first image', $
                    fid=fid1, dims=dims1, pos=pos1, /mask, m_fid=m_fid
if (fid1 eq -1) then  begin
   print,'Cancelled'
   return
endif
envi_file_query, fid1, fname=fname1, xstart=xstart1, ystart=ystart1, interleave=interleave1, $
                 nb=nb, ns=ns, nl=nl, sname=sname
print,'first image: ',fname1

if n_elements(pos1) lt 3 then begin
   print, 'Aborted'
   Message, 'Spectral subset size must be at least 3'
endif

if interleave1 eq 0 then begin
    answer = dialog_message(file_basename(fname1)+' will be coverted to BIP. Continue?',/question)
    if answer eq 'No' then begin
       print,'Cancelled'
       return
    endif
    if strmid(sname,0,3) eq '[Me' then begin
       void = dialog_message('Conversion not possible, data in memory',/error)
       print,'Conversion to BIP not possible, cancelled
       return
    endif
    dims =  [-1L,0,ns-1,0,nl-1]
    pos = lindgen(nb)
    ENVI_DOIT, 'CONVERT_INPLACE_DOIT', fid=fid1, o_interleave=2, dims=dims, pos=pos, r_fid=r_fid
    fid1 = r_fid
    interleave1 = 2
endif

; tie point
map_info = envi_get_map_info(fid=fid1)
envi_convert_file_coordinates, fid1, dims1[1], dims1[3], e, n, /to_map
map_info.mc= [0D,0D,e,n]

envi_select, title='Choose second image', fid=fid2, dims=dims2,pos=pos2
if (fid2 eq -1) then begin
   print,'Cancelled'
   return
endif
envi_file_query, fid2, fname=fname2, xstart=xstart2, ystart=ystart2, interleave=interleave2, $
                 nb=nb, ns=ns, nl=nl, sname=sname
print,'second image: ',fname2

if (n_elements(pos1) ne n_elements(pos2)) or $
   ((dims1[2]-dims1[1]) ne (dims2[2]-dims2[1])) or $
   ((dims1[4]-dims1[3]) ne (dims2[4]-dims2[3])) then begin
   print, 'Spectral/spatial subset sizes are different. Aborting.'
   Message, 'Spectral/spatial subsets different'
endif

if interleave2 eq 0 then begin
    answer = dialog_message(file_basename(fname2)+' will be coverted to BIP. Continue?',/question)
    if answer eq 'No' then begin
       print,'cancelled'
       return
    endif
    if strmid(sname,0,3) eq '[Me' then begin
       void = dialog_message('Conversion not possible, data in memory',/error)
       print,'Conversion to BIP not possible, cancelled
       return
    endif
    dims =  [-1L,0,ns-1,0,nl-1]
    pos = lindgen(nb)
    ENVI_DOIT, 'CONVERT_INPLACE_DOIT', fid=fid2, o_interleave=2, dims=dims, pos=pos, r_fid=r_fid
    fid2 = r_fid
    interleave2 = 2
endif

num_cols = dims1[2]-dims1[1]+1
num_rows = dims1[4]-dims1[3]+1
num_pixels = (num_cols*num_rows)
num_bands = n_elements(pos1)
if nbands ne num_bands then message, 'Spectral subset mismatch'

ones = fltarr(num_cols,1)+1.0  
means1 = ones##means1
means2 = ones##means2

chi_sqr=fltarr(num_cols)

; MAD output
base = widget_auto_base(title='MAD Output')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result_mad = auto_wid_mng(base)
if not result_mad.accept then $
  print, 'MAD output cancelled'

; CV output
base = widget_auto_base(title='CV Output (root name)')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result_cv = auto_wid_mng(base)
if not result_cv.accept then print, 'CV output cancelled'

print, 'Canonical correlations:'
print, reverse(rho)

vMs = reverse(sigma^2)      
sigMads=(1+fltarr(1,num_cols))##reverse(sigma)

; output to file or memory
txt =''
if result_mad.accept then begin
   if result_mad.outf.in_memory then begin
      MAD_array = fltarr(num_cols,num_rows,num_bands+1)
      txt = 'MADs -> memory,'
   end else begin
      openw, unit_mad, result_mad.outf.name, /get_lun
      txt = 'MADs -> file,'
   endelse
endif
if result_cv.accept then begin
   if result_cv.outf.in_memory then begin
      CV1_array = fltarr(num_cols,num_rows,num_bands)
      CV2_array = fltarr(num_cols,num_rows,num_bands)
      txt = txt+' CVs -> memory'
   end else begin
      openw, unit_cv1, result_cv.outf.name+'_1', /get_lun
      openw, unit_cv2, result_cv.outf.name+'_2', /get_lun
      txt = txt+' CVs -> file'
   endelse
endif
txt=txt+' ...'

tile_id1 = envi_init_tile(fid1,pos1,num_tiles=num_tiles, $
    interleave=interleave1,xs=dims1[1],xe=dims1[2],ys=dims1[3],ye=dims1[4])
tile_id2 = envi_init_tile(fid2,pos2, $
    interleave=interleave2,xs=dims2[1],xe=dims2[2],ys=dims2[3],ye=dims2[4])
progressbar = Obj_New('progressbar', Color='blue', Text='0',$
                       title=txt,xsize=250,ysize=20)
progressbar->start
abort=0
for tile_index=0L,num_tiles-1 do begin
   tile1 = envi_get_tile(tile_id1,tile_index)
   tile2 = envi_get_tile(tile_id2,tile_index)
   if interleave1 eq 1 then tile1 = transpose(tile1)
   if interleave2 eq 1 then tile2 = transpose(tile2)
   CV1s = (tile1-means1)##A
   CV2s = (tile2-means2)##B
   MADs = reverse(CV1s - CV2s,1)
   chi_sqr = float(total((MADs/sigMADs)^2,1))
   if result_mad.accept then begin
      if result_mad.outf.in_memory then begin
         for i=0,num_bands-1 do MAD_array[*,tile_index,i] = float(MADs[i,*])
         MAD_array[*,tile_index,num_bands] = chi_sqr
      end else begin
        writeu,unit_mad,[[float(MADs)],float(transpose(chi_sqr))]
      endelse
   endif
   if result_cv.accept then begin
      if result_cv.outf.in_memory then begin
         for i=0,num_bands-1 do begin
            CV1_array[*,tile_index,i] = float(CV1s[i,*])
            CV2_array[*,tile_index,i] = float(CV2s[i,*])
         endfor
      end else begin
         writeu,unit_cv1,float(CV1s)
         writeu,unit_cv2,float(CV2s)
      endelse
   endif
   if progressbar->CheckCancel() then begin
      print,'output aborted'
      abort=1
      tile_index=num_tiles
   endif
   pct=tile_index*100/num_tiles
   progressbar->Update,pct,text=strtrim(pct,2)+'%'
endfor
progressbar->destroy
envi_tile_done, tile_id1
envi_tile_done, tile_id2
if result_mad.accept then if not result_mad.outf.in_memory then free_lun, unit_mad
if result_cv.accept then if not result_cv.outf.in_memory then begin
   free_lun, unit_cv1
   free_lun, unit_cv2
endif
if abort then return

if result_mad.accept then begin
   bnames = strarr(num_bands+1)
   bnames[0:num_bands-1] = 'MAD('+file_basename(fname1)+':'+file_basename(fname2)+') '+strtrim(lindgen(num_bands)+1,2)
   bnames[num_bands] = 'CHI_SQR('+file_basename(fname1)+':'+file_basename(fname2)+')'
   if result_mad.outf.in_memory then begin
      envi_enter_data, MAD_array, $
        r_fid=r_fid, $
        map_info=map_info, $
        bnames=bnames, $
        xstart=xstart1+dims1[1], ystart=ystart1+dims1[3], $
        descrip='MAD: t1= '+file_basename(fname1)+' t2= '+file_basename(fname2)
      print, 'MADs written to memory'
      envi_assign_header_value, fid=r_fid, keyword='varMADs',value=vMs
   end else begin
      envi_setup_head ,fname=result_mad.outf.name, ns=num_cols, $
        r_fid=r_fid, $
        nl=num_rows, nb=num_bands+1, $
        data_type=4, interleave=2, /write, /open, $
        bnames=bnames, $
        map_info=map_info, $
        xstart=xstart1+dims1[1], ystart=ystart1+dims1[3], $
        descrip='MAD: t1= '+fname1+' t2= '+fname2
      print, 'File created ', result_mad.outf.name
      envi_assign_header_value, fid=r_fid, keyword='varMADs',value=vMs
      envi_write_file_header,r_fid
   endelse
endif

envi_assign_header_value, fid=r_fid, keyword='varMADs', $
    value=2*(1-rho)

if result_cv.accept then begin
   bnames1 = strarr(num_bands)
   bnames2 = strarr(num_bands)
   bnames1= 'CV1('+file_basename(fname1)+':'+file_basename(fname2)+') '+strtrim(lindgen(num_bands)+1,2)
   bnames2= 'CV2('+file_basename(fname1)+':'+file_basename(fname2)+') '+strtrim(lindgen(num_bands)+1,2)
   if result_cv.outf.in_memory then begin
      envi_enter_data, CV1_array, $
        map_info=map_info, $
        bnames=bnames1, $
        xstart=xstart1+dims1[1], ystart=ystart1+dims1[3], $
        descrip='CV1: t1= '+file_basename(fname1)+' t2= '+file_basename(fname2)
      envi_enter_data, CV2_array, $
        map_info=map_info, $
        bnames=bnames2, $
        xstart=xstart1+dims1[1], ystart=ystart1+dims1[3], $
        descrip='CV2: t1= '+file_basename(fname1)+' t2= '+file_basename(fname2)
      print, 'CVs written to memory'
   end else begin
      envi_setup_head ,fname=result_cv.outf.name+'_1', ns=num_cols, $
        nl=num_rows, nb=num_bands, $
        data_type=4, interleave=2, /write, /open, $
        bnames=bnames1, $
        map_info=map_info, $
        xstart=xstart1+dims1[1], ystart=ystart1+dims1[3], $
        descrip='CV1: t1= '+fname1+' t2= '+fname2
      envi_setup_head ,fname=result_cv.outf.name+'_2', ns=num_cols, $
        nl=num_rows, nb=num_bands, $
        data_type=4, interleave=2, /write, /open, $
        bnames=bnames2, $
        map_info=map_info, $
        xstart=xstart1+dims1[1], ystart=ystart1+dims1[3], $
        descrip='CV2: t1= '+fname1+' t2= '+fname2
      print, 'File created ', result_cv.outf.name+'_1'
      print, 'File created ', result_cv.outf.name+'_2'
   endelse
endif

End