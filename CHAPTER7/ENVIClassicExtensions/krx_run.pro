; docformat = 'rst'
; krx_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO krx_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Kernel RX (CUDA)', $
      REF_VALUE = 'SAM Target Finder with BandMax', $
      EVENT_PRO = 'krx_run', $
      UVALUE = 'KRX',$
      POSITION = 'after' 
END

function eps, X, double=double
  if n_elements(double) eq 0 then double = 0
  ma = machar(double=double)
  return, X*float((ma.ibeta))^(ma.machep)
end

;+
; :Description:
;      Performs kernel RX anomaly detection
;      using a Gaussian kernel::      
;         Kwon et al (2005). 
;         Kernel RX-algorithm: A nonlinear
;         anomaly detector for hyperspectral
;         imagery. IEEE TGARS 43(2) 388-397      
; :Params:
;      event:  in, required 
;         if called from the ENVI menu                  
; :Uses:
;      ENVI, CENTER, KERNEL_MATRIX, GPUKERNEL_MATRIX, 
;      COYOTE, GPULIB, CLUST_WTS1
; :Author:
;      Mort Canty (2013)       
;-
pro KRX_run, event 

COMPILE_OPT IDL2

seed = 12345L
  
print, '-------------------'
print, 'Kernel RX'
print, systime(0)
print, '-------------------'

catch, theError
if theError ne 0 then begin
   void = Dialog_Message(!Error_State.Msg, /error)
   return
endif

envi_select, title='Choose hyper- or multispectral image', $
             fid=fid, dims=dims,pos=pos
if (fid eq -1) then begin
   print, 'cancelled'
   return
endif
envi_file_query, fid, fname=fname,xstart=xstart,ystart=ystart
num_cols = dims[2]-dims[1]+1
num_rows = dims[4]-dims[3]+1
num_bands = n_elements(pos)
num_pixels = num_cols*num_rows
print, 'input file '+fname

; sample size
base = widget_auto_base(title='Sample size, 0 for k-means determination')
wg = widget_sslider(base, title='Samples', min=0, max=5000L, $
  value=1000L, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif 
m = long(result.slide)
if m gt 0 then $
   print,'Training on '+strtrim(m,2)+' samples' $
else print, 'Using k-means ...'   

; nscale parameter for Gaussian kernel
base = widget_auto_base(title='NSCALE')  
we = widget_param(base, dt=4, field=3, floor=0.,xsize= 50,$  
  default=10.0, uvalue='param', /auto)  
result = auto_wid_mng(base)  
if (result.accept ne 0) then nscale = float(result.param) $
   else nscale = 10.0 
   
; output destination
base = widget_auto_base(title='KRX Output')
sb = widget_base(base, /row, /frame)
wp = widget_outfm(sb, uvalue='outf', /auto)
result1 = auto_wid_mng(base)
if (result1.accept eq 0) then begin
  print, 'Output cancelled'
  return
end

; image data matrix
GG = dblarr(num_bands,num_pixels)
for i=0,num_bands-1 do $
 GG[i,*] = envi_get_data(fid=fid,dims=dims,pos=pos[i]) 

if m gt 0 then begin
; random training data matrix  
   indices = randomu(seed,m,/long) mod num_pixels
   G = GG[*,indices]
endif else begin
; k-means training matrix
   m = 100
   start_time = systime(2)    
   G = CLUST_WTS1(GG,N_CLUSTERS=m)
   print,'Elapsed time for k-means: '+strtrim(systime(2)-start_time,2)+' seconds'   
endelse  

start_time = systime(2)  

K = kernel_matrix(G,gma=gma,nscale=nscale)
Kc = center(K)
print,'GMA: ' + strtrim(gma,2)

print, 'Pseudoinverse of centered kernel matrix ...'
; pseudoinvert centered kernel matrix
lambda = la_eigenql(Kc,/double,eigenvectors=alpha)
alpha = transpose(alpha)
idx = reverse(sort(lambda))
lambda = lambda[idx]
alpha = alpha[idx,*]
tol = m*eps(max(lambda))
_ = where(lambda gt tol, r)
alpha = alpha[0:r,*]
lambda = lambda[0:r]
Kci = alpha##diag_matrix(1.0/lambda)##transpose(alpha)
print, 'Nonzero single precision eigenvalues: '+strtrim(r,2)
print,'Elapsed time: '+strtrim(systime(2)-start_time,2)+' seconds'

start_time = systime(2)  

print,'Calculating anomaly image ...'
progressbar = Obj_New('cgprogressbar',$
             title='Kernel RX ',xsize=300,ysize=20,/cancel)            
progressbar->start
if gpu_detect() then begin
  outimage_gpu = gpuFltArr(num_cols,num_rows)
  GG_gpu = gpuputarr(float(GG))   
  G_gpu = gpuputarr(float(G))
  K_gpu = gpuputarr(float(K)) 
  Kci_gpu = gpuputarr(float(Kci))
  onesm_gpu = gpuputarr(fltarr(m)+1.0)
  onesnc_gpu = gpuputarr(fltarr(num_cols)+1.0)  
  tmp_gpu = gpuTotal(K_gpu,1) 
  a = -1.0*gpuTotal(K_gpu)/m^2   
  tmp1_gpu = gpuAdd((1.0/m),tmp_gpu,0.0,tmp_gpu,a) 
  Ku_gpu = gpumatrix_multiply(tmp1_gpu,onesnc_gpu,/btranspose)    
  gpuFree, [tmp_gpu,tmp1_gpu]
  i = 0L
  while i lt num_rows do begin
     pct=i*100/num_rows
     if progressbar->CheckCancel() then begin
         print,'RX aborted'
         i = num_rows-1
     endif     
     progressbar->update,fix(pct)
     gpuView,GG_gpu,i*num_cols*num_bands,num_bands*num_cols,GGi_gpu 
     gpuReform,GGi_gpu,num_bands,num_cols             
     tmp_gpu = gpukernel_matrix(GGi_gpu,G_gpu,gma=gma,nscale=nscale)      
     a_gpu = gpuTotal(tmp_gpu,1)    
     a1_gpu = gpumatrix_multiply(onesm_gpu,a_gpu,/btranspose)         
     Kg_gpu = gpuAdd(1.0,tmp_gpu,(-1.0/m),a1_gpu,0.0)   
     gpuFree, [tmp_gpu,a_gpu,a1_gpu]
     Kgu_gpu = gpuSub(Kg_gpu,Ku_gpu)        
     gpuView,outimage_gpu,i*num_cols,num_cols,rowi_gpu          
     rowi_gpu = gputotal(Kgu_gpu*(Kgu_gpu##Kci_gpu),1,LHS=rowi_gpu)    
     gpuFree, [Kgu_gpu,Kg_gpu]
     i++ ; next row
   endwhile
   outimage = float(gpuGetarr(outimage_gpu))
end else begin 
   outimage = fltarr(num_cols,num_rows)
   Ku = total(K,1)/m-total(K)/m^2   
   Ku = (dblarr(num_cols)+1.0)##Ku  
   i = 0L                   
   while i lt num_rows do begin
      pct=i*100/num_rows  
      if progressbar->CheckCancel() then begin
         print,'RX aborted'
         i = num_rows-1
      endif           
      progressbar->update,fix(pct)      
      GGi = GG[*,i*num_cols:(i+1)*num_cols-1]
      Kg = kernel_matrix(GGi,G,gma=gma,nscale=nscale)     
      a = total(Kg,1,/double)                   
      a = a##(dblarr(m)+1.0)                     
      Kg = Kg - a/m                                 
      Kgu = Kg - Ku                         
      d = total(Kgu*(Kgu##Kci),1,/double)
      outimage[*,i] = float(d)
      i++ ; next row
   endwhile
endelse    

progressbar->Destroy     
print,'Elapsed time: '+strtrim(systime(2)-start_time,2)+' seconds'  
  
; map tie point
map_info = envi_get_map_info(fid=fid)
envi_convert_file_coordinates,fid,dims[1],dims[3],e,n,/to_map
map_info.mc = [0D,0D,e,n]

; write to memory or file
bnames=['anomalies']
if (result1.outf.in_memory eq 1) then begin
   envi_enter_data, outimage, $
      map_info=map_info, $
      bnames=bnames, $
      xstart=xstart+dims[1], ystart=ystart+dims[3], $
      descrip='Kernel RX: '+file_basename(fname)
   print, 'Result written to memory'
end else begin
   openw, unit, result1.outf.name, /get_lun
   writeu, unit, outimage
   envi_setup_head,fname=result1.outf.name, ns=num_cols, nl=num_rows, nb=1, $
                    data_type=4, $
                    interleave=0, $
                    file_type=0, $
                    map_info=map_info, $
                    xstart=xstart+dims[1], $
                    ystart=ystart+dims[3], $
                    bnames=bnames,$            
                    descrip='Kernel RX: '+file_basename(fname), $
                    /write,/open
   print, 'File created ', result1.outf.name
   free_lun, unit
endelse

end
