pro EX5_2
; create a LoG filter (and plot it)
   sigma = 2.0
   filter = fltarr(16,16)
   for i=0L,15 do for j=0L,15 do filter[i,j] = $
   (1/(2*!pi*sigma^6))*((i-8)^2+(j-8)^2-2*sigma^2) $
    *exp(-((i-8)^2+(j-8)^2)/(2*sigma^2))
;   thisDevice =!D.Name
;   set_plot, 'PS'
;   Device, Filename='fig5_4.eps',xsize=3,ysize=3, $
;            /encapsulated
   shade_surf,filter
;   device,/close_file
;   set_plot,thisDevice
; get an image band and pad it
   envi_select, title='Choose multispectral band', $
             fid=fid, dims=dims,pos=pos, /band_only
   num_cols = dims[2]-dims[1]+1
   num_rows = dims[4]-dims[3]+1
   image = fltarr(num_cols+16,num_rows+16)
   image[0:num_cols-1,0:num_rows-1] = $
       envi_get_data(fid=fid,dims=dims,pos=0)
; pad the filter as well
   filt = image*0
   filt[0:15,0:15] = filter
; perform filtering in frequency domain
   image= float(fft(fft(image)*fft(filt),1))
; get zero-crossing positions
   indices = where( (image*shift(image,1,0) lt 0) $
              or (image*shift(image,0,1) lt 0) )
; perform Sobel filter for edge strengths
   edge_strengths = sobel(image)
; create the contour image and return it to ENVI
   image = image*0
   image[indices]=edge_strengths[indices]
   envi_enter_data, image
end