pro ex9_2
; kernelized PCA for nonlinear change detection
; between two single-band images

m =   1000L   ; training sample size
gma = 1/5.   ; radial basis kernel parameter
seed= 12345

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

; full design matrix
GG = fltarr(2,num_pixels)
GG[0,*] = envi_get_data(fid=fid1,dims=dims1,pos=0)
GG[1,*] = envi_get_data(fid=fid2,dims=dims2,pos=0)

; artificial nonlinearity
GG[1,*] = (GG[1,*])^3

GG = standardize(GG)

; sample design matrix
indices = randomu(seed,m,/long) mod num_pixels
G = GG[*,indices]

; linear principal components transformation
C = correlate(GG,/covariance,/double)
C = (C+transpose(C))/2
void = eigenql(C, eigenvectors=W, /double)
PCs = GG##transpose(W)
image = fltarr(num_cols,num_rows,2)
for i=0,1 do image[*,*,i]=reform(PCs[i,*],num_cols,num_rows)
envi_enter_data, image

; radial basis kernel Gram matrix
K = fltarr(m,m)
for i=0,m-1 do for j=i,m-1 do begin
   ex = G[*,i]-G[*,j]
   K[i,j]=(K[j,i]=exp(-gma*total(ex^2)))
endfor

; eigenvalues and eigenvectors
lambda = eigenql(K, eigenvectors=V,/double)
print,lambda[0:9]

; dual variables
alpha = fltarr(m,m)
for i=0,m-1 do alpha[*,i] = V[*,i]/sqrt(lambda[i])

; first 10 nonlinear principal components, row-by-row
PCs = fltarr(num_cols,num_rows,10)
i = 0L
progressbar = Obj_New('progressbar', Color='blue', Text=' ',$
              title='projecting...',xsize=300,ysize=20)
progressbar->start
while i lt num_rows do begin
   pct=i*100/num_rows
   progressbar->Update,fix(pct),text=strtrim(pct,2)+'%'
; get the ith row
   g1 = GG[*,i*num_cols:(i+1)*num_cols-1]
; inner products
   exG   = exp(-gma*total(G^2,1))
   exg1  = exp(-gma*total(g1^2,1))
   exGg1 = exp(2*gma*G##transpose(g1))
   k = (transpose(exG)##exg1)*exGg1
; project
   PCs[*,i,*] = alpha[*,0:9]##k
   if progressbar->CheckCancel() then begin
      print,'interrupted...'
      i=num_rows-1
   endif
   i++
endwhile
progressbar->destroy
envi_enter_data, PCs

shade = fltarr(400,400,10)
ind = round((G+3)*400./6.) > 0
ind = transpose(ind) < 399
g1 = fltarr(2,400)
g1[0,*] = findgen(400)*6./400.-3.
for i=0,399 do begin
   g1[1,*] = i*6./400.-3.
   shade[*,i,0:1] = W##transpose(g1)
endfor
black = min(shade)
for i=0,1 do begin
  shadei = shade[*,*,i]
  shadei[ind[*,0],ind[*,1]]=black
  shade[*,*,i]=shadei
endfor
envi_enter_data, reverse(shade[*,*,0:1],2)

for i=0,399 do begin
   g1[1,*] = i*6./400.-3.
   exG   = exp(-gma*total(G^2,1))
   exg1  = exp(-gma*total(g1^2,1))
   exGg1 = exp(2*gma*G##transpose(g1))
   k = (transpose(exG)##exg1)*exGg1
   shade[*,i,*] = alpha[*,0:9]##k
endfor
black = min(shade)
for i=0,9 do begin
  shadei = shade[*,*,i]
  shadei[ind[*,0],ind[*,1]]=black
  shade[*,*,i]=shadei
endfor
envi_enter_data, reverse(shade,2)

end
