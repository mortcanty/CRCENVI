pro ex8_1
; extended K-means

envi_select, title='Choose multispectral band', $
             fid=fid, dims=dims,pos=pos, /band_only
if (fid eq -1) then return
cols = dims[2]-dims[1]+1
rows = dims[4]-dims[3]+1
m = float(cols*rows)

; map tie point
map_info = envi_get_map_info(fid=fid)
envi_convert_file_coordinates, fid, $
   dims[1], dims[3], e, n, /to_map
map_info.mc[2:3]= [e,n]

; image band
G = bytscl(envi_get_data(fid=fid,dims=dims,pos=0))
; classification image
labels = G*0

; parameters
sigma2 = stddev(float(G)-shift(G,[1,0]))^2
alphaE = -1/(2*alog(1.0/8))

; initial K, means and priors 
hist = histogram(G,nbins=256)
indices = where(hist,K)
means = (findgen(256))[indices]
priors = hist[indices]/m

; iteration
progressbar = Obj_New('progressbar', Color='blue', $
  title='EKM clustering: delta...',xsize=250,ysize=20)
progressbar->start
delta = 100.0 & iter = 0
while (delta gt 1.0) and (iter le 100) do begin
   if progressbar->CheckCancel() then begin
      print,'clustering aborted'
      progressbar->Destroy & return
   endif
   progressbar->Update,iter,text=strtrim(delta,2)
; drop (almost) empty clusters
   indices = where(priors gt 0.01, K)
; distance array   
   ds = fltarr(m,K)  
   means = means[indices]
   priors = priors[indices]
   means1 = means
   priors1 = priors
   means = means*0.0
   priors = priors*0.0
   for j=0,K-1 do begin
      ms = replicate(means1[j],m)
      logps = replicate(alog(priors1[j]),m)
      ds[*,j] = (G-ms)^2/(2*sigma2) - alphaE*logps
   endfor
   min_ds = min(ds,dimension=2)
   for j=0,K-1 do begin
      indices = where(ds[*,j] eq min_ds,count)
      if count gt 0 then begin
         mj = n_elements(indices)
         priors[j] = mj/m
         means[j] = total(G[indices])/mj
         labels[indices] = j+1
      endif
   endfor
   delta = max(abs(means - means1))
   iter++
endwhile
progressbar->Destroy

envi_enter_data, labels, file_type=3, $
  num_classes=K+1, $
  class_names='class '+strtrim(indgen(K+1),2), $
  lookup=class_lookup_table(indgen(K+1)), $
  map_info=map_info

end
