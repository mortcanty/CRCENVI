pro ex2_1
   thisDevice =!D.Name
   set_plot, 'PS'
   Device,Filename='fig2_1.eps',xsize=6,ysize=4,$
      /inches,/encapsulated
   plot, histogram(total(randomu(seed,10000,12),2), $
          nbins=12)
   device,/close_file
   set_plot, thisDevice  
end