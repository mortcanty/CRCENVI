pro ex3_4
envi_select, title='Choose multispectral image', $
             fid=fid, dims=dims,pos=pos
if (fid eq -1) then begin
   print, 'Canceled'
   return
endif
if n_elements(pos) lt 2 then begin
   print, 'Aborted'
   Message, 'Spectral subset size must be at least 2'
endif
envi_file_query, fid, fname=fname
cols = dims[2]-dims[1]+1
rows = dims[4]-dims[3]+1
bands = n_elements(pos)
; data matrix for difference image
D = fltarr(bands,cols*rows)
for i=0,bands-1 do begin
   temp=float(envi_get_data(fid=fid,dims=dims, $
              pos=pos[i]))
   D[i,*]=temp-(shift(temp,1,0)+shift(temp,0,1))/2.0
endfor
; noise covariance
S_N = correlate(D,/covariance,/double)   
print,'Noise covariance matrix, file ', $
     file_basename(fname)
print,S_N
end