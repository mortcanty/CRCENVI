; docformat = 'rst'
; maxlike_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details. 
     
PRO maxlike_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'MaxLike (test output)', $
      REF_VALUE = 'Maximum Likelihood', $
      EVENT_PRO = 'maxlike_run', $
      UVALUE = 'MAXLIKE',$
      POSITION = 'after'
END    
    
;+
; :Description:
;       ENVI extension for classification of
;       a multispectral image with maximum likelihood
;       with option for reserving a random sample
;       of training data for subsequent evaluation    
; :Params:
;      event:  in, required 
;         if called from the ENVI menu                    
; :Uses:
;      ENVI::
;      DIFFERENCE
; :Author:
;      Mort Canty (2009)       
;-
Pro maxlike_run, event

COMPILE_OPT IDL2

catch, theError
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
   return
endif


   print, '---------------------------------'
   print, 'Maximum Likelihood Supervised Classification'
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
   base = widget_auto_base(title='Output MaxLike classification to file')
   sb = widget_base(base, /row, /frame)
   wp = widget_outf(sb, uvalue='outf', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
     print, 'output of classification aborted'
     return
   endif
   class_filename=result.outf
   openw, class_unit, class_filename, /get_lun
   base = widget_auto_base(title='Output MaxLike rule image to file')
   sb = widget_base(base, /row, /frame)
   wp = widget_outf(sb, uvalue='outf', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
      print, 'output of rule image cancelled'
      rule_filename='none'
      rule_flag=0
   end else begin
      rule_filename=result.outf
      rule_flag=1
      openw, rule_unit, rule_filename, /get_lun
   endelse

; output file for test results

  outfile=dialog_pickfile(filter='*.tst',/write,/overwrite_prompt,title='Save test results to disk')
  if (outfile eq '') then begin
     print,'cancelled'
     return
  endif

; construct the training examples 
   Gs = fltarr(num_bands)
   Ls = fltarr(K)
   for i=0,K-1 do begin
      Gs1 = envi_get_roi_data(roi_ids[ptr[i]],$
          fid=fid, pos=pos)
      Gs  = [[Gs], [Gs1]]
      Ls1 = fltarr(K,n_elements(Gs1[0,*]))
      Ls1[i,*]=1.0
      Ls = [[Ls], [Ls1]]
   endfor
   Gs = Gs[*,1:n_elements(Gs[0,*])-1]
   Ls = Ls[*,1:n_elements(Ls[0,*])-1]
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
   
; train the classifier
   mn =  fltarr(num_bands,K)
   cov = fltarr(num_bands,num_bands,K)
   class_names = strarr(K+1)
   class_names[0]='unclassified'
   void = max(transpose(Ls),labels,dimension=2)
   labels = byte((labels/m))
   for i=0,K-1 do begin
      class_names[i+1]='class'+string(i+1)
      indices = where(labels eq i,count)
      if count gt 1 then begin
         GGs = Gs[*,indices]
         for j=0,num_bands-1 do mn[j,i] = mean(GGs[j,*])
         cov[*,*,i] = correlate(GGs,/covariance)
      endif
   endfor

envi_check_save, /classification

   lookup = [[0,0,0],[roi_colors[*,ptr]]]
   if rule_flag then begin
       envi_doit, 'class_doit', fid=fid, dims=dims, $
            pos=pos, out_bname='MaxLike('+fname+')', $
            class_names=class_names, $
            lookup=lookup, $
            method=2, mean=mn, cov=cov, $
            out_name = class_filename, $
            rule_out_name = rule_filename, $
            rule_out_bname=roi_names[ptr]
       print, 'classification file created ', class_filename
       print, 'rule file created ', rule_filename
   end  else  begin
      envi_doit, 'class_doit', fid=fid, dims=dims, $
            pos=pos, out_bname='MaxLike('+fname+')', $
            class_names=class_names, $
            lookup=lookup, $
            method=2, mean=mn, cov=cov, $
            out_name = class_filename
      print, 'classification file created ', class_filename
   endelse

; test the classifier
   dummy = fltarr(1,num_test,num_bands)
   dummy[0,*,*] = transpose(Gs_test)
   dummy_dims = [-1L,0,0,0,num_test-1]
   dummy_pos = indgen(num_bands)
   envi_enter_data, dummy, r_fid=r_fid
   envi_doit, 'class_doit', fid=r_fid, r_fid=r_fid1, $
         dims=dummy_dims, pos=dummy_pos, out_bname='Test', $
         class_names=class_names, lookup=lookup, $
         method=2, mean=mn, cov=cov, /in_memory
   test_classes = (envi_get_data(fid=r_fid1,pos=0,dims=dummy_dims))[*]
   envi_file_mng, id=r_fid,/remove
   envi_file_mng, id=r_fid1,/remove
   openw,lun,outfile,/get_lun
   printf,lun,'; MaxLike test results for '+ fname
   printf,lun,'; '+systime(0)
   printf,lun,'; Classification image: '+ class_filename
   printf,lun,'; Class rule image: '+ rule_filename
   printf,lun, num_test, K
   void = max(transpose(Ls_test),labels,dimension=2)
   labels = byte((labels/num_test)+1)
   printf,lun, transpose([[test_classes], [labels]])
   free_lun,lun
   print, 'test results written to ' + outfile
  
End