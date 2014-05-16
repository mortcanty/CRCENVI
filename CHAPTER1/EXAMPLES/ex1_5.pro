pro EX1_5
   COMPILE_OPT IDL2
   
   envi_select, title='Choose multispectral image', $
                fid=fid, dims=dims,pos=pos
   if (fid eq -1) then begin
      print, 'cancelled'
      return
   endif
   envi_file_query, fid, fname=fname
   cols = dims[2]-dims[1]+1
   rows = dims[4]-dims[3]+1
   bands = n_elements(pos)
   pixels = cols*rows
; data matrix
   G = fltarr(bands,pixels)
; read in image and subtract means
   for i=0,bands-1 do begin
     temp = envi_get_data(fid=fid,dims=dims,pos=pos[i])
     G[i,*] = temp-mean(temp)
   end
; principal components transformation
   C = correlate(G,/covariance,/double)
   void = eigenql(C, eigenvectors=U, /double)
   PCs = G##transpose(U)
; BSQ array
   image = fltarr(cols,rows,bands)
   for i = 0, bands-1 do image[*,*,i] = $
      reform(PCs[i,*],cols,rows)
; map tie point
   map_info = envi_get_map_info(fid=fid)
   envi_convert_file_coordinates, fid, $
      dims[1], dims[3], e, n, /to_map
   map_info.mc = [0D,0D,e,n]
   envi_enter_data, image, map_info = map_info
end

