function chirp, t, t0
  result = t-t
  idx = where((t-t0 lt 2000) and (t-t0 gt 0)) 
  tt = t[idx]-t0
  result[idx] = sin(2*!pi*2e-3*(tt+1e-3*tt^2))
  return, result
end

pro EX4_2
  thisDevice = !D.Name
  set_plot, 'PS'
  Device, FileName = 'fig4_2.eps',xsize=5,ysize=4, $
          /inches,/encapsulated   
  t = findgen(5000L)
  plot, t, chirp(t,400)+12, yrange = [-8.0,12.0],$
                           xtitle = 'Time'
  oplot, chirp(t,800)+8
  oplot, chirp(t,1400)+4
  signal = chirp(t,400)+chirp(t,800)+chirp(t,1400)
  kernel = (chirp(t,0))[0:2000]  
  oplot, signal
  oplot, 0.003*convol(signal,kernel)-6
  Device, /close_file
  set_plot,thisDevice
end

