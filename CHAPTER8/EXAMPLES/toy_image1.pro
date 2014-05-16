PRO toy_image1
 image = fltarr (400 ,400 ,2)
 n = randomu(seed ,400 ,400,/normal )
 n1 = 8*randomu(seed, 400, 400)-4
 image[*,*,0] = n1+8
 image[*,*,1]=n1^2+0.3*randomu(seed,400,400,/normal )+8
 image[0:199,*,0]=randomu(seed, 200,400,/normal )/2+8
 image[0:199,*,1]=randomu(seed, 200,400,/normal )+14
 envi_enter_data,image
end