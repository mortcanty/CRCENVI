; docformat = 'rst'
; kpca_run5.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

pro kpca_run5_extensions_init
   compile_opt idl2
   e5 = envi(/current)
   e5.addextension, 'Kernel PCA (CUDA)', 'kpca_run5', path='Chapter4'
end

;+
; :Description:
;      Performs kernel principal components analysis
;      using a Gaussian kernel::      
;         Nielsen, A. A. and Canty, M. J. (2008). 
;         Kernel principal component analysis
;         for change detections. In SPIE Europe 
;         Remote Sensing Conference, Cardiff,
;         Great Britain, 15-18 September, volume 7109.       
; :Params:
;      None                  
; :Uses:
;      ENVI, CENTER, KERNEL_MATRIX, GPUKERNEL_MATRIX, GPUKERNELPROJECT, COYOTE, GPULIB
; :Author:
;      Mort Canty (2013)          
;-
pro KPCA_run5

COMPILE_OPT idl2

catch, err
   if err ne 0 then begin
      catch, /cancel
      if e5 ne !null then $
         e5.reporterror, 'Error: ' + !error_state.msg $
      else $
         message, !error_state.msg, /continue, /noname
      message, /reset
      return
   endif
   
e5 = envi(/current)
   if e5 eq !null then $
      message, 'This extension requires an interactive ENVI session 
      
seed = 12345L
  
print, '-------------------'
print, 'Kernel PCA'
print, systime(0)
print, '-------------------'

UI = e5.UI
inRaster = UI.SelectInputData(/raster,/disable_sub_rect,bands=bands,title='Choose image')
; get image bands in BIP format
if inRaster ne !Null then data = inRaster.GetData(interleave='bip') else return         

sz = size(data)
num_bands = sz[1]
num_cols = sz[2]
num_rows = sz[3]
num_pixels = num_cols*num_rows

; sample size
base = widget_auto_base(title='Sample size, 0 for k-means')
wg = widget_sslider(base, title='Samples', min=0, max=5000, $
  value=2000, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif 
m = long(result.slide)
if m gt 0 then $
   print,'Training on '+strtrim(m,2)+' samples' $
else print, 'Using k-means'    

; number of PCs
base = widget_auto_base(title='Numer of PCs')
wg = widget_sslider(base, title='PCs', min=5, max=50, $
  value=10, dt=1, uvalue='slide', /auto)
result = auto_wid_mng(base)
if (result.accept eq 0) then begin
   print, 'cancelled'
   return
endif 
num_pcs = long(result.slide)

; nscale parameter for Gaussian kernel
base = widget_auto_base(title='NSCALE')  
we = widget_param(base, dt=4, field=3, floor=0.,xsize= 50,$  
  default=1.0, uvalue='param', /auto)  
result = auto_wid_mng(base)  
if (result.accept ne 0) then nscale = float(result.param) $
   else nscale = 1.0 
   
query=dialog_message('Center on training means',/question)
if query eq 'Yes' then ct=1 else ct=0   

progressbar = Cgprogressbar(/nocancel,title='Initializing...')               
progressbar->start

; image data matrix
GG = reform(data,num_bands,num_pixels)

if m gt 0 then begin
; random training data matrix  
   indices = randomu(seed,m,/long) mod num_pixels
   G = GG[*,indices]
endif else begin
; k-means training matrix
   m = 100
   G = CLUST_WTS(GG,N_CLUSTERS=m)
   print,'Elapsed time for k-means: '+strtrim(systime(2)-start_time,2)+' seconds'   
endelse  

start_time = systime(2) 

; centered radial basis kernel matrix
if gpu_detect() then begin
   print,'running CUDA ...'
   G_gpu = gpuputarr(G)
   K_gpu = gpukernel_matrix(G_gpu,gma=gma,nscale=nscale)
   K = gpugetarr(K_gpu)
   gpufree,[G_gpu,K_gpu]
end else begin
   print,'CUDA not available ...'
   K = kernel_matrix(G,gma=gma,nscale=nscale)
endelse     
K = center(K)

print,'GMA: ' + strtrim(gma,2)

progressbar->Update,0

; eigenvalues and eigenvectors of centered kernel matrix
print, 'diagonalizing kernel...'
lambda = la_eigenql(K,/double,eigenvectors=V,range=[m-num_pcs,m-1])
idx = reverse(sort(lambda))
lambda = lambda[idx]
V = V[*,idx]

print,'First '+strtrim(num_pcs,2)+' eigenvalues of the kernel matrix'
print, lambda
print, 'Fraction of variance explained in first 3 PCS'
print, total(lambda[0:2])/total(lambda)
   
envi_plot_data,findgen(num_pcs)+1,lambda[0:num_pcs-1],$
   plot_title='Kernel PCA',$
   title='Kernel PCA',$
   ytitle='Eigenvalue',$
   xoff=100,yoff=100

print,'Elapsed time for training: '+strtrim(systime(2)-start_time,2)+' seconds'
progressbar->destroy

start_time = systime(2) 

; dual variables (normalized eigenvectors)
alpha = float(diag_matrix(1/sqrt(lambda))##V) 
   
if not gpuKernelProject(alpha,  gma, G, GG, num_cols, num_rows, image, center_train=ct) then $
       message, 'projection aborted'

print,'Elapsed time for projection: '+strtrim(systime(2)-start_time,2)+' seconds'

; save the PCs as a raster
tempFile = e5.GetTemporaryFilename()
outRaster = e5.CreateRaster(tempFile,image, $
  inherits_from=inRaster,interleave='bsq')
outRaster.Save

print, 'Kernel PCs saved to temporary file.'

view = e5.GetView()
layer1 = view.CreateLayer(inRaster,bands=[0,1,2])
layer2 = view.CreateLayer(outRaster,bands=[0,1,2])
       
end
