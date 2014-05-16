pro HSV5
COMPILE_OPT IDL2
e = envi()
; get MS image
envi_select, title='Select 3-band MS input file', $
                    fid=fid1,dims=dims1,pos=pos1
if (fid1 eq -1) or (n_elements(pos1) ne 3) then return
; get PAN image
envi_select, title='Select panchromatic image', $
             fid=fid2,pos=pos2,dims=dims2,/band_only
if (fid2 eq -1) then return
; linear stretch the images and convert to byte format
envi_doit,'stretch_doit',fid=fid1,dims=dims1,pos=pos1,$
  method=1, r_fid=r_fid1,out_min=0,out_max=255, $
  range_by=0,i_min=0,i_max=100,out_dt=1,/in_memory
envi_doit,'stretch_doit',fid=fid2,dims=dims2,pos=pos2,$
  method=1,r_fid=r_fid2,out_min=0,out_max=255, $
  range_by=0,i_min=0,i_max=100,out_dt=1,/in_memory
envi_file_query,r_fid2,ns=f_ns,nl=f_nl
f_dims = [-1l,0,f_ns-1,0,f_nl-1]
; HSV pan-sharpening
envi_doit, 'sharpen_doit', $
  fid=[r_fid1,r_fid1,r_fid1],pos=[0,1,2],f_fid=r_fid2,$
  f_dims=f_dims,f_pos=[0],method=0,interp=0,$
  r_fid=r_fid3,/in_memory
out_raster = ENVIFIDToRaster(r_fid3)
view = e.GetView()
layer = view.CreateLayer(out_raster)

end