; docformat = 'rst'
; svm_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
    

PRO svm_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'SVM (test output)', $
      REF_VALUE = 'Support Vector Machine', $
      EVENT_PRO = 'svm_run', $
      UVALUE = 'SVM',$
      POSITION = 'after'
END    

;+
; :Description:
;       ENVI extension for running and testing
;       classification of a multispectral image with SVM.
;       Kernel parameters can only be set programmatically. 
;       Classification and rule images are written to memory.             
; :Uses:
;      ENVI
;      DIFFERENCE
; :Params:
;      event:  in, required 
;         if called from the ENVI menu ;      
; :Author:
;      Mort Canty (2013)      
;-
Pro svm_run, event

COMPILE_OPT IDL2

catch, theError
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
   return
endif

   print, '---------------------------------'
   print, 'SVM Supervised Classification'
   print, systime(0)
   print, '---------------------------------'
; select the image to be classified
   envi_select, title='Enter File for classification', fid=fid, pos=pos, dims=dims
   if (fid eq -1) then begin
      print,'cancelled'
      return
   endif
   envi_file_query, fid, fname=fname, xstart=xstart, ystart=ystart
   print, 'file: '+fname
   num_cols=dims[2]-dims[1]+1
   num_rows=dims[4]-dims[3]+1
   n_bands=n_elements(pos)

; SVM parameters (assumes radial basis (2) kernel type only):
   kernel_type = 2 
   pyramid_levels = 2
   kernel_gamma = 1.0/n_bands
   penalty = 100.0

; output file for test results
  outfile=dialog_pickfile(filter='*.tst',/write,/overwrite_prompt,title='Save test results to disk')
  if (outfile eq '') then begin
     print,'cancelled'
     return
  endif

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
   
   roi_addr = envi_get_roi(roi_ids[ptr[0]])
   roi_label = intarr(n_elements(roi_addr)) 
   for i=1,K-1 do begin
      tmp = envi_get_roi(roi_ids[ptr[i]])
      roi_addr  = [roi_addr,tmp]
      roi_label = [roi_label,intarr(n_elements(tmp))+i]  
   endfor
   m = n_elements(roi_addr)   
; split into training and test data
   seed = 12345L
   num_test = m/3
; sampling with replacement      
   test_indices = randomu(seed,num_test,/long) mod m     
   test_indices =  test_indices[sort(test_indices)]
   train_indices=difference(lindgen(m),test_indices)
   m = n_elements(train_indices)
   roi_addr_train = roi_addr[train_indices]
   roi_addr_test  = roi_addr[test_indices]  
   roi_label_train = roi_label[train_indices]
   roi_label_test  = roi_label[test_indices]  
   
; delete ROIs   
   envi_save_rois,'roi_save$$',envi_get_roi_ids()
   envi_get_roi_information, roi_ids, roi_colors=roi_colors
   envi_delete_rois, envi_get_roi_ids()     
   
; replace with training ROIs
   for i=0, K-1 do begin
      r = roi_colors[0,i]
      g = roi_colors[1,i]
      b = roi_colors[2,i]
      color=envi_color_match(r,g,b)
      roi_id = envi_create_roi(ns=num_cols,nl=num_rows, $
                color=color[0], name='region '+strtrim(i,2))
      idx = where(roi_label_train eq i,count)
      if count gt 0 then begin
         pts = roi_addr_train[idx]
         envi_define_roi, roi_id,/point,xpts=pts mod num_cols,ypts=pts/num_cols
      endif   
   endfor
   roi_ids = envi_get_roi_ids(fid=fid)
   n_ids = n_elements(roi_ids)
   roi_ids = roi_ids[1:n_ids-1] ; strip off the empty ROI

; run the SVM classifier
   envi_doit, 'envi_svm_doit', $
    fid=fid, pos=pos, dims=dims, $
    roi_ids=roi_ids, $
    /in_memory,r_fid=r_fid, $
    /rule_in_memory, $
    kernel_type=kernel_type, $
    pyramid_levels=pyramid_levels, $
    kernel_gamma=kernel_gamma, $
    penalty=penalty
    
; get the classification image for test
   class_image = envi_get_data(fid=r_fid,pos=0,dims=dims)
   
; test the classifier
   openw,lun,outfile,/get_lun
   printf,lun,'; SVM test results for '+ fname
   printf,lun,'; '+systime(0)
   printf,lun,'; Classification image: not generated'
   printf,lun,'; Class rule image: not generated'
   printf,lun, num_test, K
   for i=0,K-1 do begin
      idx = where(roi_label_test eq i,count)
      if count gt 0 then begin
         idx1 = roi_addr_test[idx] 
         SVM_labels = class_image[idx1]
         labels = roi_label_test[idx]+1
         printf,lun,transpose([[SVM_labels], [labels]])    
      endif
   endfor
   free_lun, lun
   print, 'Test results written to '+outfile
   
; resore the original ROI data
   envi_delete_rois,envi_get_roi_ids()
   envi_restore_rois, 'roi_save$$'
   file_delete, 'roi_save$$', /noexpand_path

End