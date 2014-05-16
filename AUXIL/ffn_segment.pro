Pro ffn_segment
   print, '---------------------------------------------'
   print, 'FFN Segment Classification with Hu moments'
   print, systime(0)
   print, '---------------------------------------------'

   FFN='FFNKal'

; select the segment data file
   filename=dialog_pickfile(filter='*.csv',title='Select segment data')
   if filename eq '' then return
   print, 'Reading training data from '+filename
   openr,unit,filename,/get_lun
; get first line
   aLine=''
   readF,unit,aLine
   header = strsplit(aLine,';',count=count,/extract)
   num_features=count-10
   features=fltarr(num_features+1)
   while ~ EOF(unit) do begin
      readf,unit,aLine
      entries=strsplit(aLine,';',/extract)
      entries=[entries[2],entries[10:10+num_features-1]]
      features=[[features],[float(entries)]]
   endwhile
   free_lun,unit
   num_segments=(size(features))[2]-1
   features=features[*,1:num_segments]
; select segment image file
   envi_select,fid=fid,pos=pos,dims=dims,/no_dims,/no_spec,/band_only,title='Select segment image'
   if fid eq -1 then return
   segment_image = envi_get_data(fid=fid,dims=dims,pos=pos)
; select band image
   no_band=0
   envi_select,fid=fid1,pos=pos1,dims=dims1,/no_dims,/no_spec,/band_only,title='Select image band (optional)'
   if fid1 eq -1 then  no_band=1 else band_image = envi_get_data(fid=fid1,dims=dims1,pos=pos1)
   num_cols = dims[2]-dims[1]+1
   num_rows = dims[4]-dims[3]+1
; calculate table of invariant moments
   moments = fltarr(7,num_segments)
   progressbar = Obj_New('progressbar', Color='blue', Text=' ',$
                    title='Calculating moments...',xsize=300,ysize=20)
   progressbar->start
   window, 11,xsize=400,ysize=400
   A0=fltarr(5,5)
   for i=0,num_segments-1 do begin
      if progressbar->CheckCancel() then begin
         progressbar->destroy
         wdelete,11
         return
      endif
      pct = (i)*100.0/(num_segments)
      progressbar->Update,fix(pct)
      indices = where(segment_image eq i,count)
      if count gt 0 then begin
         X = indices mod num_cols
         Y = indices/num_cols
         A = segment_image[min(X):max(X),min(Y):max(Y)]
         indices1 = where(A eq i,complement=complement,ncomplement)
         if no_band then A[indices1]=1.0 else begin
            B = band_image[min(X):max(X),min(Y):max(Y)]
            A[indices1]=B[indices1]
         endelse
         if ncomplement gt 0 then A[complement]=0.0
         if i mod 10 eq 0 then begin
             wset,11
             tvscl,A0
             tvscl,A
             A0=A*0
         endif
         moments[*,i]=hu_moments(A)
      endif
   endfor
   progressbar->destroy
   wdelete,11

   features=[[transpose(features)],[transpose(moments)]]

;   features=[[transpose(features[0,*])],[transpose(moments)]]


; build up example pairs
   indices = where(features[*,0] ne -1)
   labels = features[indices,0]
   num_classes=max(labels)
   Xs = features[indices,1:(size(features))[2]-1]
   Ys = fltarr(n_elements(indices),num_classes)
   for i=0,num_classes-1 do begin
      indices = where(labels eq i+1,count)
      if count gt 0 then Ys[indices,i]=1.0
   endfor
   max_x = max(Xs,dimension=1)
   min_x = min(Xs,dimension=1)
   for i=0,(size(Xs))[2]-1 do Xs[*,i] = (Xs[*,i]-min_x[i])/(max_x[i]-min_x[i])
; train the network
   base = widget_auto_base(title='Number of hidden units')
   wg = widget_sslider(base, title='Hidden units', min=2, max=16, $
     value=4, dt=1, uvalue='slide', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then n_hu = 4 else  n_hu = result.slide
   print, 'Training ...'
   FFN = Obj_New(FFN,Xs,Ys,n_hu)
   FFN->train
   print, '... done'
; classify
   print, 'classifying...'
   Xs = features[*,1:(size(features))[2]-1]
   for i=0,(size(Xs))[2]-1 do Xs[*,i] = (Xs[*,i]-min_x[i])/(max_x[i]-min_x[i])
   classes = FFN->classify(Xs,Probs)
   Obj_Destroy, FFN
   progressbar = Obj_New('progressbar', Color='blue', Text=' ',$
                    title='Classifying segments...',xsize=300,ysize=20)
   progressbar->start
   for i=0,num_segments do begin
      if progressbar->CheckCancel() then begin
         progressbar->destroy
         return
      endif
      pct = (i)*100.0/(num_segments)
      progressbar->Update,fix(pct)
      indices= where(segment_image eq i,count)
      if count gt 0 then segment_image[indices]=classes[i]
   endfor
   envi_enter_data,segment_image,file_type=envi_file_type('ENVI Classification'), $
      class_names=string(indgen(num_classes+1)), $
      num_classes=num_classes+1, $
      lookup=class_lookup_table(indgen(num_classes+1))
   progressbar->destroy
   print, '...done.'

End