pro EX3_1
envi_select, title='Choose multispectral band', $
             fid=fid, dims=dims,pos=pos, /band_only
if (fid eq -1) then begin
   print, 'cancelled'
   return
endif
cols = dims[2]-dims[1]+1
rows = dims[4]-dims[3]+1
image = envi_get_data(fid=fid,dims=dims,pos=pos[0])
; arrays of i and j values
a = lindgen(cols,rows)
i = a mod cols
j = a/cols
; shift Fourier transform to center
image = (-1)^(i+j)*image
; compute power spectrum an return to ENVI
envi_enter_data, alog((abs(FFT(image)))^2)
end

