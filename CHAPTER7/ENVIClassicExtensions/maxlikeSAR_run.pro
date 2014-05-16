; docformat = 'rst'
; maxlikeSAR_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details. 
     
PRO maxlikeSAR_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'MaxLike polSAR', $
      REF_VALUE = 'Maximum Likelihood', $
      EVENT_PRO = 'maxlikeSAR_run', $
      UVALUE = 'C_CORRECTION',$
      POSITION = 'after'
END    
    
;+
; :Description:
;       ENVI extension for classification of
;       a polarimetric SAR image with maximum likelihood
;       with option for reserving a random sample
;       of training data for subsequent evaluation    
; :Params:
;      event:  in, required 
;         if called from the ENVI menu                    
; :Uses:
;      ENVI
; :Author:
;      Mort Canty (2014)       
;-
Pro maxlikeSAR_run, event

COMPILE_OPT STRICTARR

   print, '---------------------------------------------------'
   print, 'Maximum Likelihood polSAR Supervised Classification'
   print, systime(0)
   print, '---------------------------------------------------'
      
;   catch, theError
;   if theError ne 0 then begin
;      void = Dialog_Message(!Error_State.Msg, /error)
;      return
;   endif   
   
; input multi-look averaged covariance image
   envi_select, title='Choose (spatial subset of) covariance matrix image', $
     fid=fid, dims=dims, pos=pos
   if (fid eq -1) then begin
     print, 'cancelled'
     return
   end
   envi_file_query, fid, fname=fname, bnames=bnames,xstart=xstart, ystart=ystart
   cols = dims[2]-dims[1]+1
   rows = dims[4]-dims[3]+1
   bands = n_elements(pos)
   ; bands = 9: quad pol
   ; bands = 4  dual pol
   ; bands = 1  single pol
   
; tie point
   map_info = envi_get_map_info(fid=fid)
   envi_convert_file_coordinates, fid, dims[1], dims[3], e, n, /to_map
   map_info.mc= [0D,0D,e,n]

; read in image as data matrix
   GG1s = fltarr(bands,cols*rows)
   for i = 0,bands-1 do GG1s[i,*] = envi_get_data(dims=dims, fid=fid, pos=pos[i])   
; get associated ROIs
   envi_delete_rois, /all
   fname1 = envi_pickfile(filter='*.roi',title='Select ROI data file')
   envi_restore_rois, fname1
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

; output file for test results
  testfile=dialog_pickfile(filter='*.tst',/write,/overwrite_prompt,title='Save test results to disk')
  if (testfile eq '') then begin
     print,'cancelled'
     return
  endif  
  
; make complex data matrix
  case bands of
    9: begin
      c_bands = 6
      GGs = complexarr(6,cols*rows)
      GGs[0,*] = complex(GG1s[0,*],0.0*GG1s[0,*])
      GGs[1,*] = complex(GG1s[1,*],GG1s[2,*])
      GGs[2,*] = complex(GG1s[3,*],GG1s[4,*])
      GGs[3,*] = complex(GG1s[5,*],0.0*GG1s[5,*])
      GGs[4,*] = complex(GG1s[6,*],GG1s[7,*])
      GGs[5,*] = complex(GG1s[8,*],0.0*GG1s[8,*])
    end
    4: begin
      c_bands = 3
      GGs = complexarr(3,m)
      GGs[0,*] = complex(GG1s[0,*],0.0*GG1s[0,*])
      GGs[1,*] = complex(GG1s[1,*],GG1s[2,*])
      GGs[2,*] = complex(GG1s[3,*],0.0*GG1s[3,*])
    end
  endcase  

; construct the training examples 
   Gs = fltarr(c_bands)
   Ls = fltarr(K)
   for i=0,K-1 do begin    
      idx = envi_get_roi(roi_ids[ptr[i]])
      Gs1 = fltarr(c_bands,n_elements(idx))
      for j=0,c_bands-1 do Gs1[j,*] = GGs[pos[j],idx]      
      Gs  = [[Gs], [Gs1]]
      Ls1 = fltarr(K,n_elements(Gs1[0,*]))
      Ls1[i,*]=1.0
      Ls = [[Ls], [Ls1]]
   endfor
   envi_file_mng, /remove, id=fid    
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
   Gst = Gs[*,test_indices]
   Lst = Ls[*,test_indices]
   m = n_elements(train_indices)
   Gs = Gs[*,train_indices]
   Ls = Ls[*,train_indices]
   
; train the classifier
   void = max(transpose(Ls),labels,dimension=2)
   labels = byte((labels/m))
   D = fltarr(K) 
   case bands of
   9: begin  
;       quad polarimetric    
         Sinv = complexarr(9,K)
         for i=0,K-1 do begin
            indices = where(labels eq i,count)
            if count gt 1 then begin
               G = mean(Gs[*,indices],dimension=2)          
               Sinv[*,i] = invert( [[G[0],      G[1],      G[2]],$
                                   [conj(G[1]), G[3],      G[4]],$
                                   [conj(G[2]),conj(G[4]),G[5]]] )
               K1    = G[0]
               A1    = G[1]
               RHO1  = G[2]
               XSI1  = G[3]
               B1    = G[4]
               ZETA1 = G[5]
               D[i] = k1*xsi1*zeta1 + 2*real_part(a1*b1*conj(rho1)) - xsi1*(abs(rho1)^2) - k1*(abs(b1)^2) - zeta1*(abs(a1)^2) 
            endif
         endfor
      end
   4: begin     
;       dual polarimetric
         Sinv = complexarr(4,K)
         for i=0,K-1 do begin
            indices = where(labels eq i,count)
            if count gt 1 then begin
               G = mean(Gs[*,indices],dimension=2)     
               Sinv[*,i] = invert( [[G[0],       G[1]],$
                                    [conj(G[1]), G[2]]] )
               K1    = G[0]
               A1    = G[1]
               XSI1  = G[2]
               D[i] =  k1*xsi1 - abs(a1)^2 
            endif
         endfor
       end 
   endcase            
; test the classifier
   discr = fltarr(K,num_test)
   case bands of
      9: tmp = [ Gst[0,*], conj(Gst[1,*]), conj(Gst[2,*]), Gst[1,*], Gst[3,*], conj(Gst[4,*]), Gst[2,*], Gst[4,*], Gst[5,*] ]
      4: tmp = [ Gst[0,*], conj(Gst[1,*]), Gst[1,*], Gst[2,*] ]
   endcase   
   for i = 0,K-1 do begin
      SSi = transpose(fltarr(num_test)+1)##Sinv[*,i]
      discr[i,*] = -alog(D[i]) - total(SSi*tmp,1) 
   endfor
   void = max(transpose(discr),test_classes,dimension=2)
   test_classes = byte((test_classes/num_test)+1)
   openw,lun,testfile,/get_lun
   printf,lun,'; MaxLike test results for '+ fname
   printf,lun,'; '+systime(0)   
   printf,lun,'; Classification image: N/A'
   printf,lun,'; Class probilities image: N/A'   
   printf,lun, num_test, K
   void = max(transpose(Lst),labels,dimension=2)
   labels = byte((labels/num_test)+1)
   printf,lun, transpose([[test_classes], [labels]])
   free_lun,lun
   print, 'test results written to ' + testfile
; classify the image    
   discr = fltarr(K,cols*rows)
   case bands of
      9: tmp = [ GGs[0,*], conj(GGs[1,*]), conj(GGs[2,*]), GGs[1,*], GGs[3,*], conj(GGs[4,*]), GGs[2,*], GGs[4,*], GGs[5,*] ]
      4: tmp = [ GGs[0,*], conj(GGs[1,*]), GGs[1,*], GGs[3,*] ]
   endcase   
   for i = 0,K-1 do begin
      SSi = transpose(fltarr(cols*rows)+1)##Sinv[*,i]
      discr[i,*] = -alog(D[i]) - total(SSi*tmp,1) 
   endfor
   void = max(transpose(discr),classes,dimension=2)   
   classes = reform(byte(classes/(cols*rows)+1),cols,rows)
   lookup = [[0,0,0],[roi_colors[*,ptr]]]   
   envi_enter_data, classes, $
                    data_type=1, $
                    file_type=3, $
                    map_info=map_info,$
                    num_classes=K+1, $
                    class_names=['unclassified',roi_names[ptr]], $
                    lookup = lookup, $
                    bnames=['maxlikeSAR('+fname+')']
    
End