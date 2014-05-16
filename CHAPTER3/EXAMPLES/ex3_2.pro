function psi_m, x
   if x lt 0.0 then return, 0.0
   if x lt 0.5 then return, 1.0
   if x lt 1.0 then return, -1.0
   return, 0.0
end
function psi, m, k, n
   c = 2^n
   result = fltarr(c)
   x = findgen(c)/c
   for i=0,c-1 do result[i]=psi_m(2^m*x[i]-k)
   return, result
end

pro ex3_2
; generate wavelet basis B_8
   n = 8L
   B = fltarr(2^n,2^n)+1.0
   i = 1
   for m=0,n-1 do for k=0,2^m-1 do begin
      B[i,*] = psi(m,k,n)
      i++
   endfor
; get a 256x256 grayscale image
   envi_select, title='Choose multispectral band', $
                fid=fid, dims=dims,pos=pos, /band_only
   G = envi_get_data(fid=fid,dims=dims,pos=pos)
   G = float(G[0:255,0:255])
   print, 'Size original image:',256L*256L*4L,' bytes'
   envi_enter_data, G + 0.0
; transform the columns and rows
   for i=0,255 do G[i,*] = invert(B)##G[i,*]
   for j=0,255 do G[*,j] = invert(B)##transpose(G[*,j])
   envi_enter_data, G + 0.0
; convert to sparse format
   G = sprsin(G,thresh=1.0)
   write_spr, G, 'sparse.dat'
   openr, 1,'sparse.dat' & status = fstat(1) & close,1
   print, 'Size compressed image:',status.size,' bytes'
; invert the transformation
   G = fulstr(G)
   for j=0,255 do G[*,j] = B##transpose(G[*,j])
   for i=0,255 do G[i,*]  = B##G[i,*]
   envi_enter_data, G
end

