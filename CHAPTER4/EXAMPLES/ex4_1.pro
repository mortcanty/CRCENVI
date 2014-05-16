pro EX4_1

; get an image band from ENVI
   envi_select, title='Choose multispectral band', $
             fid=fid, dims=dims,pos=pos, /band_only
   if (fid eq -1) then return
   num_cols = (c = dims[2]-dims[1]+1)
   num_rows = dims[4]-dims[3]+1
; pick out the center row of pixels
   image = envi_get_data(fid=fid,dims=dims,pos=pos[0])
   g = float(image[*,num_rows/2])
; define a FIR kernel of length m = 5
   h = [1,2,3,2,1]
; setup for postscript 
   thisDevice = !D.Name
   set_plot, 'PS'
   Device, FileName = 'fig4_1.eps',xsize=5,ysize=3, $
          /inches,/encapsulated             
; convolve in the spatial domain and plot       
   plot, convol(g,h,center=0)
; pad the arrays to c + m - 1
   g = [g,[0,0,0,0]]
   hp = g*0
   hp[0:4] = h
; convolve in the frequency domain and plot
   oplot, fft(c*fft(g)*fft(hp),1)-200
   Device, /close_file
   set_plot,thisDevice

end

