pro subset_copy
   COMPILE_OPT IDL2  
   
   dim = 256L
   
; get an image band for T1
   envi_select, title='Choose T1 image', $
                fid=fid1, dims=dims1,pos=pos1
   if (fid1 eq -1) then return
   
; get an image band for T2
   envi_select, title='Choose T2 image', $
             fid=fid2, dims=dims2,pos=pos2
   if (fid2 eq -1) then return
   
; read images
   num_cols = dims1[2]-dims1[1]+1
   num_rows = dims1[4]-dims1[3]+1
   num_bands = n_elements(pos1)
   num_pixels = num_cols*num_rows
   im1 = fltarr(num_cols,num_rows,num_bands)
   im2 = im1*0.0
   for i=0,num_bands-1 do begin
     im1[*,*,i] = envi_get_data(fid=fid1,dims=dims1,$
                  pos=pos1[i])
     im2[*,*,i] = envi_get_data(fid=fid2,dims=dims2, $
                  pos=pos2[i])
   endfor
   
; copy subset with Gaussian noise
   im2[0:dim-1,0:dim-1,*] = im1[0:dim-1,0:dim-1,*] $
      + 0.1*randomu(seed,dim,dim,num_bands,/normal)
   envi_enter_data, im2

end
