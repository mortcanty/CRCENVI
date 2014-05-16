; docformat = 'rst'
; hcl_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO hcl_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Agglomerative Hierachical', $
      REF_VALUE = 'K-Means', $
      EVENT_PRO = 'hcl_run', $
      UVALUE = 'HCL',$
      POSITION = 'after'
END

;+
; :Description:
;       ENVI extension for agglomerative hierarchical
;       clustering
; :Params:
;       event:  in, optional 
;          required if called from ENVI                 
; :Uses:
;       ENVI::
;       HCL::  
;       CLASS_LOOKUP_TABLE::            
; :Author:
;       Mort Canty (2013)      
;-
pro HCL_run, event

COMPILE_OPT IDL2

print, '---------------------------------'
print, 'Agglomerative Hierarchic Clustering'
print, systime(0)
print, '---------------------------------'

catch, theError
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
   return
endif

envi_select, title='Choose multispectral image for clustering', $
             fid=fid, dims=dims,pos=pos
if (fid eq -1) then begin
   print, 'cancelled'
   return
endif

envi_file_query, fid, fname=fname
print, 'Selected image: ',fname
print, 'Selected bands: ',pos + 1

; number of samples
base = widget_auto_base(title='Number of Samples')
wg = widget_sslider(base, title='Samples', min=500, max=2000, $
  value=1000, dt=2, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif
n_samples = result.slide

; number of clusters
base = widget_auto_base(title='Number of Classes')
wg = widget_sslider(base, title='Classes', min=2, max=16, $
  value=6, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif
K = byte(result.slide)
print, 'Number of Classes', K

; output destination
base = widget_auto_base(title='HCL Output')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
  print, 'cancelled'
  return
endif

num_cols = dims[2]-dims[1]+1
num_rows = dims[4]-dims[3]+1
num_pixels = (num_cols*num_rows)
num_bands = n_elements(pos)

; randomly sampled data matrix
indices = randomu(seed,n_samples,/long) mod num_pixels
G = fltarr(num_bands,n_elements(indices))
for i=0,num_bands-1 do begin
    temp= envi_get_data(fid=fid,dims=dims,pos=pos[i])
    G[i,*]=temp[indices]
endfor

; run hcl
print, 'Hierarchic clustering algorithm running on ',n_elements(indices),' pixels ...'
HCL, G, K, Ls

answer=dialog_message('Do you want to use present results to classify a larger spatial subset?', /question)
if answer eq 'Yes' then begin
   envi_select, title='Choose a new spatial subset', fid=fid1, dims=dims1, /no_spec
   if (fid1 ne fid) then begin
      print,'Different file. Continuing with original subset.'
      Message, 'Different file. Continuing with original subset.'
   end else dims=dims1
endif

; supervised classification
print,'Supervised classification with clustered samples ...'

mn =  fltarr(num_bands,K)
cov = fltarr(num_bands,num_bands,K)
class_names = strarr(K+1)
class_names[0]='unclassified'
for i=0,K-1 do begin
   class_names[i+1]='cluster'+string(i+1)
   indices = where(Ls eq i,count)
   if count gt 1 then begin
      for j=0,num_bands-1 do mn[j,i] = mean(G[j,indices])
      cov[*,*,i] = correlate(G[*,indices],/covariance)
   endif
endfor

envi_check_save, /classification

if (result.outf.in_memory eq 1) then in_memory=1 else in_memory=0
envi_doit, 'class_doit', fid=fid, dims=dims, pos=pos, out_bname='HCL('+fname+')', $
            class_names=class_names, lookup=class_lookup_table(indgen(K+1)), $
            method=2, mean=mn, cov=cov, in_memory=in_memory, out_name = result.outf.name
if in_memory eq 0 then print, 'File created ', result.outf.name

end


