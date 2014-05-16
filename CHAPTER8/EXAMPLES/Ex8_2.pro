Pro Ex8_2
; number of clusters
   K = 3  
; select segment file
   envi_select,fid=fid,pos=pos,dims=dims, $
     /no_dims,/no_spec,/band_only, $
     title='Select segment file'
   if fid eq -1 then return
   num_cols = dims[2]-dims[1]+1
   seg_im = envi_get_data(fid=fid,dims=dims,pos=pos) 
   num_segs = max(seg_im)          
; calculate table of invariant moments
   moments = fltarr(7,num_segs)
   for i=0L, num_segs-1 do begin
      indices = where(seg_im eq i+1,count)
      if count gt 0 then begin
         X = indices mod num_cols
         Y = indices/num_cols
         A = seg_im[min(X):max(X)+1,min(Y):max(Y)+1]
         indices1=where(A eq i+1, complement=comp,$
                           ncomplement=ncomp)
         A[indices1]=1.0 
         if ncomp gt 0 then A[comp]=0.0
         moments[*,i] = hu_moments(A,/log)         
      endif
   endfor
; agglomerative hierarchic clustering   
   HCL,standardize(moments),K,Ls   
; re-label segments with cluster label
   _ = histogram(seg_im,reverse_indices=r)
   for i=1L,num_segs do begin
       p = r[r[i]:r[i+1]-1] ;pixels in segment 
       seg_im[p]=Ls[i-1]+1  ;new label
   endfor
; output to memory   
   envi_enter_data, seg_im, $
      file_type=3,num_classes=K+1, $
      class_names=string(indgen(K+1)), $
      lookup=class_lookup_table(indgen(K+1))
End