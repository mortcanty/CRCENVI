pro EX2_2

envi_select, title='Choose multispectral image', $
             fid=fid, dims=dims,pos=pos
if (fid eq -1) then begin
   print, 'cancelled'
   return
endif

envi_file_query,fid,fname=fname,interleave=interleave
if interleave ne 2 then begin
   print, 'not BIP format'
   return
endif

num_cols = dims[2]-dims[1]+1
num_bands = n_elements(pos)

cpm = Obj_New("CPM",num_bands)

tile_id = envi_init_tile(fid,pos,num_tiles=num_tiles, $
       interleave=2,xs=dims[1],xe=dims[2], $
       ys=dims[3],ye=dims[4])

; spectral tiling
for tile_index = 0L, num_tiles-1 do $
   cpm->update, envi_get_tile(tile_id,tile_index)

print, 'Covariance matrix for image '+fname
print, cpm->covariance()

Obj_Destroy, cpm
envi_tile_done, tile_id

end

