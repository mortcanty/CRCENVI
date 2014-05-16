pro EX5_1
; create filter
   g = fltarr(512,512)
   g[0:2,0:2] = [[1,0,-1],[2,0,-2],[1,0,-1]]
; shift Fourier transform to center
   a = lindgen(512,512)
   i = a mod 512
   j = a/512
   g = (-1)^(i+j)*g
; output Fourier power spectrum as EPS file
   thisDevice =!D.Name
   set_plot, 'PS'
   Device, Filename='fig5_1.eps',xsize=3,ysize=3, $
            /inches,/encapsulated
   tvscl, (abs(FFT(g)))^2
   device,/close_file
   set_plot,thisDevice
end

