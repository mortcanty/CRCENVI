pro ex5_5

; read a JPEG image, cut out two 512x512 pixel arrays
filename = Dialog_Pickfile(Filter='*.jpg',/Read)

if filename eq '' then return else begin
   Read_JPeG,filename,image
   g1 = image[0,10:521,10:521]
   g2 = image[0,0:511,0:511]
   
; perform Fourier transforms
   f1 = fft(g1, /double)
   f2 = fft(g2, /double)
   
; Determine the offset
   g = fft(f2*conj(f1)/abs(f1*conj(f1)),$
        /inverse,/double)
   pos = where(g eq max(g))
   print,'Offset='+strtrim(pos mod 512) $
        +strtrim(pos/512)
        
; output as EPS file
   thisDevice =!D.Name
   set_plot, 'PS'
   Device, Filename='c:\temp\fig5_25.eps', $
          xsize=4,ysize=4,/inches,/Encapsulated
   shade_surf,g[0,0:50,0:50]
   device,/close_file
   set_plot, thisDevice
endelse

end