pro ex9_1

; get an image band for T1
envi_select, title='Choose multispectral band for T1',$
             fid=fid1, dims=dims1,pos=pos1, /band_only
if (fid1 eq -1) then return

; get an image band for T2
envi_select, title='Choose multispectral band for T2',$
             fid=fid2, dims=dims2,pos=pos2, /band_only
if (fid2 eq -1) then return

if (dims1[2]-dims1[1] ne dims2[2]-dims2[1]) or $
   (dims1[4]-dims1[3] ne dims2[4]-dims2[3]) then return
num_cols = dims1[2]-dims1[1]+1
num_rows = dims1[4]-dims1[3]+1
num_pixels = num_cols*num_rows
bitemp = fltarr(2,num_pixels)
bitemp[0,*] = envi_get_data(fid=fid1,dims=dims1,pos=0)
bitemp[1,*] = envi_get_data(fid=fid2,dims=dims2,pos=0)
bitemp[0,*] = bitemp[0,*] - mean(bitemp[0,*])
bitemp[1,*] = bitemp[1,*] - mean(bitemp[1,*])

; initial principal components transformation
C = correlate(bitemp,/covariance)
eivs = eigenql(C, eigenvectors=V, /double)
PCs = bitemp##transpose(V)
window,11,xsize=600,ysize=400, $
   title='Iterated principal components'
PLOT, [-1,1],[-abs(V[0,1]/V[0,0]),abs(V[0,1]/V[0,0])],$
  xrange=[-1,1],yrange=[-1,1],color=0, $
         background='FFFFFF'XL
         
; iterated principal components
covpm = Obj_New('CPM',2)
iter = 0L
niter = 5
; begin iteration
repeat begin
; determine weights via clustering
   sigma1 = sqrt(eivs[1])
   U = randomu(seed,num_pixels,2)
   unfrozen = where( abs(PCs[1,*]) ge sigma1, $
      complement = frozen)
   U[frozen,0] = 1.0
   U[frozen,1] = 0.0
   for j=0,1 do U[*,j]=U[*,j]/total(U,2)
   EM, PCs[1,*], U, Ms, Ps, Fs, T0=0, $
      unfrozen=unfrozen
; weighted covariance matrix
   covpm->update, bitemp, weights=U[*,0]
   C = covpm->covariance()
   eivs = eigenql(C, eigenvectors=V, /double)
   OPLOT, [-1,1],[-abs(V[0,1]/V[0,0]), $
       abs(V[0,1]/V[0,0])],linestyle=2, color=0
   PCs = bitemp##transpose(V)
   iter++
end until iter ge niter

; return result to ENVI
out_array = fltarr(num_cols,num_rows,2)
for i=0,1 do out_array[*,*,i] = $
    reform(U[*,i],num_cols,num_rows)
envi_enter_data, out_array

end
