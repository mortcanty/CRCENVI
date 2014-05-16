pro ex7_1

; get an image band
envi_select, title='Choose multispectral band', $
             fid=fid, dims=dims,pos=pos, /band_only
if (fid eq -1) then return
num_cols = dims[2]-dims[1]+1
num_rows = dims[4]-dims[3]+1
num_pixels = float(num_cols*num_rows)

; map tie point
map_info = envi_get_map_info(fid=fid)
envi_convert_file_coordinates, fid, $
   dims[1], dims[3], e, n, /to_map
map_info.mc[2:3]= [e,n]

; linear stretch the image and convert to byte format
envi_doit,'stretch_doit', fid=fid, dims=dims, pos=pos,$
  method=1, r_fid=r_fid, out_min=0, out_max=255, $
  range_by=0, i_min=2, i_max=98, out_dt=1, /in_memory
envi_file_query, r_fid, dims=dims
image = bytscl(envi_get_data(fid=r_fid,dims=dims,pos=0))
envi_file_mng, id=r_fid, /remove

; set up classification image
class_image = image*0

; parameters
sigma2 = stddev(float(image)-shift(image,[1,0]))^2
alphaE = -1/(2*alog(1.0/4))

; initial K, priors and means
hist = histogram(image)
indices = where(hist,K)
means = (findgen(256))[indices]
priors = hist[indices]/num_pixels

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
   progressbar->Update,(iter*100)/100, $
                  text=strtrim(delta,2)
; drop (almost) empty clusters
   indices = where(priors gt 0.0, K)
   means = means[indices]
   priors = priors[indices]
   ds = fltarr(num_pixels,K)
   means1 = means
   priors1 = priors
   means = means*0.0
   priors = priors*0.0
   for j=0,K-1 do begin
      ms = replicate(means1[j],num_pixels)
      logps = replicate(alog(priors1[j]),num_pixels)
      ds[*,j] = (image-ms)^2/(2*sigma2) - alphaE*logps
   endfor
   min_ds = min(ds,dimension=2)
   for j=0,K-1 do begin
      indices = where(ds[*,j] eq min_ds,count)
      if count gt 0 then begin
         nj = n_elements(indices)
         priors[j] = nj/num_pixels
         means[j] = total(image[indices])/nj
         class_image[indices] = j+1
      endif
   endfor
   delta = max(abs(means - means1))
   iter++
endwhile
progressbar->Destroy

envi_enter_data, class_image, file_type=3, $
  num_classes=K+1, $
  class_names='class '+strtrim(indgen(K+1),2), $
  lookup=class_lookup_table(indgen(K+1)), $
  map_info=map_info

end
