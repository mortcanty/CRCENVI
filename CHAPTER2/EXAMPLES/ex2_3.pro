function MI, image1, image2
; returns the mutual information of
; two grayscale byte images
   p12 = hist_2d(image1,image2,min1=0,max1=255, $
                 min2=0,max2=255)
   p12 = float(p12)/total(p12)
   p1  = histogram(image1,min=0,max=255)
   p1  = float(p1)/total(p1)
   p2  = histogram(image2,min=0,max=255)
   p2  = float(p2)/total(p2)
   p1p2 = transpose(p2)##p1
   i = where(p1p2 gt 0 and p12 gt 0)
   return, total(p12[i]*alog(p12[i]/p1p2[i]))
end

pro ex2_3
  envi_select, title='Choose first single band image',$
             fid=fid1, dims=dims,pos=pos1, /band_only
  if (fid1 eq -1) then return
  envi_select,title='Choose second single band image',$
             fid=fid2, dims=dims,pos=pos2, /band_only
  if (fid2 eq -1) then return
  image1 = $
   (envi_get_data(fid=fid1,dims=dims,pos=pos1))[*]
  image2 = $
   (envi_get_data(fid=fid2,dims=dims,pos=pos2))[*]
  print, MI(bytscl(image1),bytscl(image2))
end