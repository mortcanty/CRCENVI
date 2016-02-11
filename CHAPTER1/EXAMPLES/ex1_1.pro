pro EX1_1

   envi_select, title='Choose multispectral image', $
             fid=fid, dims=dims,pos=pos
   if (fid eq -1) then begin
      print, 'cancelled'
      return
   endif

   envi_file_query, fid, fname=fname

; image dimensions
   cols = dims[2]-dims[1]+1
   rows = dims[4]-dims[3]+1
   bands = n_elements(pos)

; BSQ array
   image = fltarr(cols,rows,bands)

   for i=0,bands-1 do image[*,*,i] = $
     envi_get_data(fid=fid,dims=dims,pos=pos[i])

; display first band
   window, 11, xsize=cols, ysize=rows, title=fname
   tvscl, image[*,*,0]

end

