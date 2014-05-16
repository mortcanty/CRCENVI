pro EX1_3

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

; array with pixel vectors as rows
G = fltarr(bands,pixels)

for i=0,bands-1 do $
   G[i,*]= envi_get_data(fid=fid,dims=dims,pos=pos[i])
C = correlate(G,/covariance,/double)

print, 'Covariance matrix for image '+fname
print, C

end