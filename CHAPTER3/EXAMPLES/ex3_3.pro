function F, x, i
   common refinement, c0,c1,c2,c3,c4
   if (i eq 0) then if (x eq 0) then $
      return, 1.0 else return, 0.0 else $
      return, c0*f(2*x,i-1)+c1*f(2*x-1,i-1)+ $
        c2*f(2*x-2,i-1)+c3*f(2*x-3,i-1)+c4*f(2*x-4,i-1)
end

pro ex3_3

   common refinement, c0,c1,c2,c3,c4
; refinement coefficients for Haar scaling function
   c0=1 & c1=1 & c2=0 & c3=0 & c4=0
; refinement coefficients for D4 scaling function
;  c0=(1+sqrt(3))/4 & c1=(3+sqrt(3))/4
;  c2=(3-sqrt(3))/4 & c3=(1-sqrt(3))/4 & c4=0

; fourth order approximation
   n=4
   x = findgen(4*2^n)
   ff=fltarr(4*2^n)
   for i=0,4*2^n-1 do ff[i]=F(x[i]/2^n,n)

; output as EPS file
   thisDevice =!D.Name
   set_plot, 'PS'
   Device, Filename='fig3_8.eps',xsize=3,ysize=2, $
            /inches,/encapsulated
   plot,x/2^n,ff,yrange=[-1,2] 
   device,/close_file
   set_plot,thisDevice

end