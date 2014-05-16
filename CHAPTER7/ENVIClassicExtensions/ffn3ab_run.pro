; docformat = 'rst'
; ffn3ab_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO ffn3ab_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Neural Net (Adaboost)', $
      REF_VALUE = 'Neural Net', $
      EVENT_PRO = 'ffn3ab_run', $
      UVALUE = 'FFN3AB',$
      POSITION = 'after'
END

Function seq_classify, M, Xs, ps
; weighted classification of observations Xs with M FFN classifiers
; returns labels, places weighted class membership probabilities in ps
common classifiers, FFNs,alpha,K
   nx = n_elements(Xs[*,0])
   _ = FFNs[0]->classify(Xs,ps1)
   ps = alpha[0]*ps1
   for i=1,M-1 do begin
      _ = FFNs[i]->classify(Xs,ps1)
      ps = ps + alpha[i]*ps1
   endfor
   den = total(ps,2)
   for i=0,K-1 do ps[*,i]=ps[*,i]/den
   _ = max(ps,labels,dimension=2)
   return, byte(labels/nx+1)
End

;+
; :Description:
;       ENVI extension for classification of a
;       multispectral image  with a 3-layer feed forward
;       neural network trained with Kalman filter
;       and using adaptive boosting AdaBoost.M1::
;         Canty, M. J. (2009). Boosting a fast neural 
;         network for supervised land cover classification. 
;         Computers and Geosciences, accepted.     
; :Params:
;       event:  in, required 
;          if called from ENVI                   
; :Uses:
;       ENVI::
;       DIFFERENCE::
;       FFN3KAL::
;       COYOTE
; :Author:
;       Mort Canty (2013)    
;-
Pro Ffn3AB_Run, event

COMPILE_OPT IDL2

common classifiers, FFNs,alpha,K

   catch, theError
   if theError ne 0 then begin
      print, 'Error reading file: '+fileName
      return
   endif
   
   print, '---------------------------------------------'
   print, 'FFN Kalman/AdaBoost Supervised Classification'
   print, systime(0)
   print, '---------------------------------------------'
; select the image to be classified
   envi_select, title='Enter File for classification', fid=fid, pos=pos, dims=dims
   if (fid eq -1) then begin
      print,'cancelled'
      return
   endif
   envi_file_query, fid, fname=fname,xstart=xstart, ystart=ystart, interleave=interleave, ns=ns, nl=nl, nb=nb
   if interleave eq 0 then begin
      answer = dialog_message(file_basename(fname)+' will be coverted to BIP. Continue?',/question)
      if answer eq 'No' then begin
         print,'cancelled'
         return
      endif
      dims =  [-1L,0,ns-1,0,nl-1]
      pos = lindgen(nb)
      ENVI_DOIT, 'CONVERT_INPLACE_DOIT', fid=fid, o_interleave=2, dims=dims, pos=pos, r_fid=r_fid
      fid = r_fid
      interleave = 2
   endif
   print, 'file: '+fname
   num_cols=dims[2]-dims[1]+1
   num_rows=dims[4]-dims[3]+1
   num_bands=n_elements(pos)
   print, 'bands: ', pos+1
; tie point
   map_info = envi_get_map_info(fid=fid)
   envi_convert_file_coordinates, fid, dims[1], dims[3], e, n, /to_map
   map_info.mc[2:3]= [e,n]

; get associated ROIs
   roi_ids = envi_get_roi_ids(fid=fid, roi_names=roi_names, roi_colors=roi_colors)
   if (roi_ids[0] eq -1) then begin
     error=dialog_message('No ROIs associated with the selected file',/error)
     print, 'No ROIs associated with the selected file'
     print, 'done'
     return
   endif
; compound widget for ROI selection
   base = widget_auto_base(title='ROI Selection')
   wm   = widget_multi(base, list=roi_names, uvalue='list', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
      error=dialog_message('No ROIs selected',/error)
      print, 'No ROIs selected'
      print, 'done'
      return
   endif
   ptr = where(result.list eq 1, K)

; output destination
   base = widget_auto_base(title='Output FFN classification to file')
   sb = widget_base(base, /row, /frame)
   wp = widget_outf(sb, uvalue='outf', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
     print, 'output of classification aborted'
     return
   endif
   class_filename=result.outf
   openw, class_unit, class_filename, /get_lun
   base = widget_auto_base(title='Output FFN Probabilities to file')
   sb = widget_base(base, /row, /frame)
   wp = widget_outf(sb, uvalue='outf', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
      print, 'output of probabilities cancelled'
      prob_filename='none'
      pr_flag=0
   end else begin
      prob_filename=result.outf
      pr_flag=1
      openw, prob_unit, prob_filename, /get_lun
   endelse

; output file for test results
  outfile=dialog_pickfile(filter='*.tst',/write,/overwrite_prompt,title='Save test results to disk')
  if (outfile eq '') then begin
     print,'output of test results cancelled'
     return
  endif

; construct the examples
   Gs = fltarr(num_bands)
   Ls = fltarr(K)
   for i=0,K-1 do begin
      Gs1 = envi_get_roi_data(roi_ids[ptr[i]], fid=fid, pos=pos)
      Gs  = [[Gs], [Gs1]]
      Ls1 = fltarr(K,n_elements(Gs1[0,*]))
      for j=0L,n_elements(Gs1[0,*])-1 do Ls1[i,j]=1.0
      Ls = [[Ls], [Ls1]]
   endfor
   Gs = Gs[*,1:n_elements(Gs[0,*])-1]
   Ls = Ls[*,1:n_elements(Ls[0,*])-1]
   max_x = max(Gs,dimension=2)
   min_x = min(Gs,dimension=2)
   for i=0,num_bands-1 do Gs[i,*] = 2*(Gs[i,*]-min_x[i])/(max_x[i]-min_x[i])-1
   m = n_elements(Gs[0,*])

; split into training and test data
   seed = 12345L
   num_test = m/3
; sampling with replacement      
   test_indices = randomu(seed,num_test,/long) mod m     
   test_indices =  test_indices[sort(test_indices)]
   train_indices=difference(lindgen(m),test_indices)
   Gs_test = Gs[*,test_indices]
   Ls_test = Ls[*,test_indices]
   m = n_elements(train_indices)
   Gs = Gs[*,train_indices]
   Ls = Ls[*,train_indices]
   
; get neural network parameters
   base = widget_auto_base(title='Number of epochs')
   wg = widget_sslider(base, title='epochs', min=1, max=50, $
     value=5, dt=1, uvalue='slide', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then epochs = 1 else  epochs = result.slide

   base = widget_auto_base(title='Number of hidden neurons in first layer')
   wg = widget_sslider(base, title='hidden neurons 1', min=2, max=20, $
     value=10, dt=2, uvalue='slide', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then LL1 = 10 else  LL1 = result.slide

   base = widget_auto_base(title='Number of hidden neurons in second layer')
   wg = widget_sslider(base, title='hidden neurons 2', min=2, max=20, $
     value=10, dt=2, uvalue='slide', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then LL2 = 10 else  LL2 = result.slide

   print, 'epochs: ',epochs
   print, 'hidden neurons: ',LL1,LL2

   window,11,xsize=600,ysize=400,title='Boosted FFN'
; AdaBoost.M1  ----------------------------------------
   FFNs = objarr(75)   ; sequence of FFN classifiers
   alpha = dblarr(75)  ; classifier weights
   error_train = dblarr(75)
   error_test  = dblarr(75)
   _ = max(transpose(Ls_test),labs,dimension=2)
   labels_test = byte((labs/num_test)+1)
   _ = max(transpose(Ls),labs,dimension=2)
   labels_train = byte((labs/m)+1)
   D = dblarr(m) + 1.0
   D = D/total(D)
   mm = 0 & cancel = 0
   start = systime(1)
   for mm = 0,74  do begin
      trial=0
   retry:
      print, 'Training FFN '+strtrim(mm,2)+' trial '+strtrim(trial,2)
      FFNs[mm] = Obj_New('FFN3KAL',transpose(Gs),transpose(Ls),LL1,LL2,epochs=epochs,D=D)
      FFNs[mm]->train,/verbose,cancel=cancel
      if cancel then break
      labs = FFNs[mm]->classify(transpose(Gs),_)
;    error
      h = labs ne labels_train
      indices = where(h, count)
      if count gt 0 then epsilon = total(D[indices]) else break
      print,epsilon
      beta = epsilon/(1-epsilon)
      if beta ge 1.0 then begin
         trial++
         if trial gt 5 then break
         obj_destroy,FFNs[mm]
         goto, retry
      endif
      alpha[mm] = alog(1/beta)
;    update distribution
      D = D*beta^(1-h)
      D = D/total(D)
;    get the training error thus far
      labs = seq_classify(mm+1,transpose(Gs),_)
      indices = where(labs ne labels_train, count)
      error_train[mm] = double(count)/m
;    get the test error thus far
      labs = seq_classify(mm+1,transpose(Gs_test),_)
      _ = where(labs ne labels_test, count)
      error_test[mm] = double(count)/num_test
;    display the errors
      wset,11
      x = findgen(mm+1)+1
      plot,x,error_test[0:mm],color=0,background='FFFFFF'XL,/xlog,linestyle=1,thick=2,$
            xtitle = 'number of networks', ytitle= 'misclassification rate'
      oplot,x,error_train[0:mm],color=0,linestyle=2
   endfor

; test the classifier
   test_classes = seq_classify(mm,transpose(Gs_Test),probs)
   openw,lun,outfile,/get_lun
   printf,lun,'; Boosted FFN test results for '+ fname
      printf,lun,'; '+systime(0)
      printf,lun,'; Classification image: '+ class_filename
      printf,lun,'; Class probabilities image: '+ prob_filename
      printf,lun, num_test, K
      void = max(transpose(Ls_test),labels,dimension=2)
      labels = byte((labels/num_test)+1)
      printf,lun, transpose([[test_classes], [labels]])
      free_lun,lun
      print, 'test results written to '+outfile

; classify the image
   progressbar = Obj_New('cgprogressbar', Text='0',$
                 title='Classifying tile ...',xsize=250,ysize=20)
   progressbar->start
   print,'classifying...'
; get (spectral) tile id
   tile_id = envi_init_tile(fid,pos,num_tiles=num_tiles, $
       interleave=interleave,xs=dims[1],xe=dims[2],ys=dims[3],ye=dims[4])
; start (spectral) tiling
   for tile_index=0L,num_tiles-1 do begin
      if interleave eq 2 then tile = float(transpose(envi_get_tile(tile_id,tile_index))) $
                         else tile = float(envi_get_tile(tile_id,tile_index))
      tile=float(tile)                         
      for i=0,num_bands-1 do tile[*,i] = 2*(tile[*,i]-min_x[i])/(max_x[i]-min_x[i])-1
      tile = tile>(-1.0)
      tile = tile<1.0
      if progressbar->CheckCancel() then begin
         print,'Classification aborted'
         free_lun, class_unit
         if pr_flag eq 1 then free_lun, prob_unit
         progressbar->Destroy
         for i=0,74 do Obj_Destroy, FFNs[i]
         return
      endif
      progressbar->Update,tile_index*100/num_tiles
      class_row = seq_classify(mm,tile,probs)
      writeu,class_unit,class_row
      if pr_flag eq 1 then $
 ; write probabilities in BIP, byte format
         writeu,prob_unit,bytscl(transpose(probs),min=0.0,max=1.0)
   endfor ; end  of tiling
   progressbar->destroy
   envi_tile_done, tile_id
   lookup = [[0,0,0],[roi_colors[*,ptr]]]
   envi_setup_head,fname=class_filename, ns=num_cols, $
                   nl=num_rows, nb=1, $
                   data_type=1, $
                   file_type=3, $
                   map_info=map_info,$
                   num_classes=K+1, $
                   interleave=0, /write, /open, $
                   class_names=['unclassified',roi_names[ptr]], $
                   xstart=xstart, $
                   ystart=ystart, $
                   lookup = lookup, $
                   bnames=['FFN('+fname+')']
   print, 'classification file created ', class_filename
   free_lun, class_unit
   if pr_flag eq 1 then begin
      envi_setup_head, fname=prob_filename, ns=num_cols, $
                   nl=num_rows, nb=K, $
                   map_info=map_info,$
                   data_type=1, $
                   interleave=2, $ ; BIP
                   /write, /open, $
                   xstart=xstart, $
                   ystart=ystart, $
                   bnames=roi_names[ptr]
      print, 'probabilities file created ', prob_filename
      free_lun, prob_unit
   endif

   print, 'done'
   for i=0,74 do Obj_Destroy, FFNs[i]
   print, 'elapsed time (seconds): ', systime(1)-start
   window,13
         plot,alpha,color=0,background='FFFFFF'XL,linestyle=0, $
            xtitle = 'Network', ytitle= 'Weight'
End