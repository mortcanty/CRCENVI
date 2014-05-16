pro ex5_4
COMPILE_OPT STRICTARR

height = 705.0
base = 370.0
pixel_size = 15.0
envi_select, title='Choose nadir image', $
   fid=fid1, dims=dims1, pos=pos1, /band_only
if (fid1 eq -1) then return
envi_select, title='Choose back-looking image', $
   fid=fid2, dims=dims2, pos=pos2, /band_only
if (fid2 eq -1) then return
im1 = envi_get_data(fid=fid1,dims=dims1,pos=pos1)
im2 = envi_get_data(fid=fid2,dims=dims2,pos=pos2)

n_cols = dims1[2]-dims1[1]+1
n_rows = dims1[4]-dims1[3]+1
parallax = fltarr(n_cols,n_rows)

progressbar = Obj_New('progressbar',Color='blue', $
 Text='0', title='Cross correlation, column ...', $
                  xsize=250,ysize=20)
progressbar->start
for i=7L,n_cols-8 do  begin
   if progressbar->CheckCancel() then begin
      envi_enter_data,pixel_size*parallax*(height/base)
      progressbar->Destroy
      return
   endif
   progressbar->Update,(i*100)/n_cols,text=strtrim(i,2)
   for j=25L,n_rows-26 do begin
      cim = correl_images(im1[i-7:i+7,j-7:j+7], $
             im2[i-7:i+7,j-25:j+25],xoffset_b=0, $
             yoffset_b=-25,xshift=0,yshift=25)
      corrmat_analyze,cim,xoff,yoff,maxcorr,edge, $
             plateau,magnification=2
      if (yoff lt -5) or (yoff gt 18) then $
        parallax[i,j]=parallax[i,j-1] else $
        parallax[i,j]=2*yoff
   endfor
endfor
progressbar->destroy
envi_enter_data, pixel_size*parallax*(height/base)

end
