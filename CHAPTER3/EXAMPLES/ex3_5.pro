pro ex3_5
; generate two classes
   n1 = randomu(seed,1000,/normal)
   n2 = n1 + randomu(seed,1000,/normal)
   B1 = [[n1],[n2]]
   B2 = [[n1+4],[n2]]
   image = [B1,B2]
   center_x = mean(image[*,0])
; principal components analysis
   C = correlate(transpose(image),/covariance,/double)
   void = eigenql(C, eigenvectors=U, /double)
; slopes of the principal axes
   a1 = U[1,0]/U[0,0]
   a2 = U[1,1]/U[0,1]
; scatterplot and principal axes
   thisDevice =!D.Name
   set_plot, 'PS'
   Device, Filename='D:\temp\fig3_12.eps',xsize=3, $
         ysize=3,/inches,/encapsulated
   plot, image[*,0],image[*,1],psym=4,/isotropic
   oplot,[center_x-10,center_x+10], [a1*(-10),a1*(10)]
   oplot,[center_x-10,center_x+10], [a2*(-10),a2*(10)]
   device,/close_file
   set_plot,thisDevice
end