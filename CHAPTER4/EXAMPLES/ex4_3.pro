pro EX4_3

   envi_select, title='Choose multispectral band', $
             fid=fid, dims=dims,pos=pos, /band_only
   if (fid eq -1) then begin
      print, 'cancelled' & return
   endif
   cols = dims[2]-dims[1]+1
   rows = dims[4]-dims[3]+1
; get the image band from ENVI
   g = envi_get_data(fid=fid,dims=dims,pos=pos[0])
; transform it
   g_hat = fft(g)
; create a Gauss filter in the frequency domain
   sigma = 50
   d = dist(cols,rows)
   h_hat =  exp(-d^2/sigma^2)
; output surface plot of filter as EPS file
   thisDevice =!D.Name
   set_plot, 'PS'
   Device, Filename='fig4_2.eps',xsize=3,ysize=3, $
         /inches,/encapsulated
   shade_surf, shift(h_hat,cols/2,rows/2)
   device,/close_file
   set_plot,thisDevice
; multiply, do inverse FFT and return result to ENVI
   envi_enter_data, fft(g_hat*h_hat,1)    ; low-pass
   envi_enter_data, fft(g_hat*(1-h_hat),1); high-pass

end

