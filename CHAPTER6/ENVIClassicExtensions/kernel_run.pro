; docformat = 'rst'
; kernel_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO kernel_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Gaussian Kernel (CUDA)', $
      REF_VALUE = 'Support Vector Machine', $
      EVENT_PRO = 'kernel_run', $
      UVALUE = 'KERNEL',$
      POSITION = 'after'    
END

Function Output, sigma, Hs, symm=symm 
COMPILE_OPT STRICTARR
   common examples, Gs, Gs_gpu, ells, K, m, cuda  
   if n_elements(symm) eq 0 then symm = 0
   result = fltarr(n_elements(Hs[0,*]),K)
   if cuda then begin
      if symm then Kappa_gpu = $
        gpukernel_matrix(Gs_gpu,gma=0.5/sigma^2) $
      else begin
         Hs_gpu = gpuputarr(Hs)
         Kappa_gpu = gpukernel_matrix(Gs_gpu, $
                         Hs_gpu,gma=0.5/sigma^2)
         gpufree,Hs_gpu
      endelse      
      Kappa = gpugetarr(Kappa_gpu)
      gpufree,Kappa_gpu
   end $
   else Kappa=gausskernel_matrix(Gs,Hs,gma=0.5/sigma^2)
   if symm then Kappa[indgen(m),indgen(m)] = 0.0
   for j=0,K-1 do begin
      Kpa = Kappa
      idx = where(ells ne j, ncomplement=nj)
      Kpa[*,idx] = 0
      result[*,j] = total(Kpa,2)/nj      
   endfor   
   return, result
End

Function Theta, sigma
COMPILE_OPT STRICTARR
   common examples, Gs, Gs_gpu, ells, K, m, cuda
   _ = max(output(sigma,Gs,/symm),labels,dimension=2)
   _ = where(labels/m ne ells, count)
   result = float(count)/m
   oplot,[alog10(sigma)],[result],psym=5,color=0
   return, result
End

;+
; :Description:
;       ENVI extension for classification of
;       a multispectral image with a Gassian 
;       kernel-based classifier::
;           Masters, T. (1995). Advanced Algorithms 
;           for Neural Networks, A C++ ;Sourcebook. 
;           J. Wiley and Sons
; :Params:
;      event:  in, required 
;         if called from the ENVI menu                  
; :Uses:
;       ENVI::
;       DIFFERENCE::
;       MINF_BRACKET::
;       MINF_PARABOLIC::
;       GAUSSKERNEL_MATRIX::
;       GPUKERNEL_MATRIX:: 
;       GPULIB::
;       COYOTE
; :Author:
;      Mort Canty (2013)
;-
Pro kernel_run, event

   COMPILE_OPT IDL2
     
   common examples, Gs, Gs_gpu, ells, K, m, cuda
   print, '---------------------------------'
   print, 'Kernel-based Classification'
   print, systime(0)
   print, '---------------------------------'
   
   catch, theError
   if theError ne 0 then begin
      void = Dialog_Message(!Error_State.Msg, /error)
      progressbar->destroy
      return
   endif
   
; select the image to be classified
   envi_select, title='Enter File for classification', fid=fid, pos=pos, dims=dims
   if (fid eq -1) then begin
      print,'cancelled'
      return
   endif
   envi_file_query, fid, fname=fname,xstart=xstart, ystart=ystart, interleave=interleave,ns=ns, nl=nl, nb=nb
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
   print, 'bands: ',pos+1
   num_cols=dims[2]-dims[1]+1
   num_rows=dims[4]-dims[3]+1
   num_bands=n_elements(pos)
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
   base = widget_auto_base(title='Output kernel-based classification to file')
   sb = widget_base(base, /row, /frame)
   wp = widget_outf(sb, uvalue='outf', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
     print, 'output of classification aborted'
     return
   endif
   class_filename=result.outf
   openw, class_unit, class_filename, /get_lun
   base = widget_auto_base(title='Output probability image to file')
   sb = widget_base(base, /row, /frame)
   wp = widget_outf(sb, uvalue='outf', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
      print, 'output of probability image cancelled'
      prob_filename='none'
      pr_flag=0
   end else begin
      prob_filename=result.outf
      pr_flag=1
      openw, prob_unit, prob_filename, /get_lun
   endelse

; output file for test results
  test_flag = 1
  outfile=dialog_pickfile(filter='*.tst',/write,/overwrite_prompt,title='Save test results to disk')
  if (outfile eq '') then begin
     print,'output of test results cancelled'
     test_flag = 0
  endif

; construct the training examples
   Gs = fltarr(num_bands)
   Ls = fltarr(K)
   for i=0,K-1 do begin
      Gs1 = envi_get_roi_data(roi_ids[ptr[i]], fid=fid, pos=pos)
      Gs  = [[Gs], [Gs1]]
      Ls1 = fltarr(K,n_elements(Gs1[0,*]))
      Ls1[i,*]=1.0
      Ls = [[Ls], [Ls1]]
   endfor
   Gs = Gs[*,1:n_elements(Gs[0,*])-1]
   Ls = Ls[*,1:n_elements(Ls[0,*])-1]
   m = n_elements(Gs[0,*])

   if test_flag then begin
; split into training and test data
      seed = 12345L
      num_test = m/3
;    sampling with replacement      
      test_indices = randomu(seed,num_test,/long) mod m     
      test_indices =  test_indices[sort(test_indices)]
      train_indices=difference(lindgen(m),test_indices)
      Gs_test = Gs[*,test_indices]
      Ls_test = Ls[*,test_indices]
      m = n_elements(train_indices)
      Gs = Gs[*,train_indices]
      Ls = Ls[*,train_indices]
   endif

; train the classifier

; set training class labels in common block
   _ = max(transpose(Ls),ells,dimension=2)
   ells = ells/m  
   
; set GPU array for training data if CUDA available
   cuda = 0
   if gpu_detect() then begin
      print,'Running CUDA...'
      cuda = 1
      if m gt 2000 then begin
;    too many training data 
         print, 'Subsampling to 2000 training examples' 
         idx = randomu(seed,2000,/long) mod m      
         ells = ells[idx]
         Gs = Gs[*,idx]
         m = 2000
      endif
      Gs_gpu = gpuputarr(Gs) 
   end else print,'CUDA not available ...'     

   window,12,xsize=600,ysize=400, $
     title='Bracketing, please wait ...'
   wset, 12
   widget_control, /hourglass

; bracket the minimum
   s1=0.1  & s2=10.0
   minF_bracket, s1,s2,s3, FUNC_NAME="Theta"

   window,12,xsize=600,ysize=400, $
     title='Iteration ...'
   plot,[alog10(s1),alog10(s3)],[theta(s1),theta(s3)],$
      psym=5,color=0,background='FFFFFF'XL, $
      xtitle = 'log(sigma)', $
      ytitle = 'theta'

; hunt it down
   minF_parabolic, s1,s2,s3, sigma, theta_min, $
                FUNC_NAME='Theta'

   if test_flag then begin
; test the classifier
      _ = max(output(sigma,Gs_test), test_labs, dimension=2)
      test_labs = byte(test_labs/num_test + 1)
      openw,lun,outfile,/get_lun
      printf,lun,'; Gaussian kernel test results for '+ fname
         printf,lun,'; '+systime(0)
         printf,lun,'; Classification image: '+ class_filename
         printf,lun,'; Class probabilities image: '+ prob_filename
         printf,lun, num_test, K
         void = max(transpose(Ls_test),labs,dimension=2)
         labs = byte((labs/num_test)+1)
         printf,lun, transpose([[test_labs], [labs]])
         free_lun,lun
         print, 'test results written to '+outfile
   endif

; classify the image
   progressbar = Obj_New('cgprogressbar', Text='0',/cancel,$
                 title='Classifying tile ...',xsize=250,ysize=20)
   progressbar->start
   print,'classifying...'
; get (spectral) tile id
   tile_id = envi_init_tile(fid,pos,num_tiles=num_tiles, $
       interleave=interleave,xs=dims[1],xe=dims[2],ys=dims[3],ye=dims[4])
   x = fltarr(1,num_bands)
; start (spectral) tiling
   for tile_index=0L,num_tiles-1 do begin
      if interleave eq 2 then tile = envi_get_tile(tile_id,tile_index) $
                         else tile = transpose(envi_get_tile(tile_id,tile_index))
      if progressbar->CheckCancel() then begin
         print,'Classification aborted'
         free_lun, class_unit
         if pr_flag eq 1 then free_lun, prob_unit
         progressbar->Destroy
         return
      endif
      progressbar->Update,tile_index*100/num_tiles
      PVs = output(sigma,tile)
      dens = (fltarr(1,K)+1)##total(PVs,2)
      PVs = PVs/dens
      _ = max(PVs,labs,dimension=2)
      writeu,class_unit,byte(labs/num_cols + 1)
      if pr_flag eq 1 then $
 ; write probabilities in BIP, byte format
         writeu,prob_unit,bytscl(transpose(PVs),min=0.0,max=1.0)
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
                   bnames=['KERNEL('+fname+')']
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

End