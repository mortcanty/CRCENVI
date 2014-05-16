pro EX1_6
   COMPILE_OPT IDL2
; start ENVI 5 GUI
   envi5 = ENVI()
; read image
   UI = envi5.UI
   inRaster = UI.SelectInputData(/raster, $
            /disable_sub_rect,bands=bands)  
; get image bands in BIP format   
   data = inRaster.GetData(interleave='bip', $
                                 bands=bands)
   sz = size(data)
   num_bands = sz[1]
   cols = sz[2]
   rows = sz[3]                                  
; data matrix
   G = reform(data,num_bands,cols*rows)
; subtract means   
   for i=0,num_bands-1 do G[i,*] = G[i,*] $
                            - mean(G[i,*])
; principal components transformation
   C = correlate(G,/covariance,/double)
   void = eigenql(C, eigenvectors=U, /double)
   PCs = G##transpose(U) 
   PCs = reform(PCs,num_bands,cols,rows)
; create and save the PCs as a raster   
   tempFile = envi5.GetTemporaryFilename()
   outRaster = envi5.CreateRaster(tempFile,PCs, $
            inherits_from=inRaster,interleave='bip')
   outRaster.Save
; display the original and transformed images
   view = envi5.GetView()
   layer1 = view.CreateLayer(inRaster,bands=[0,1,2])
   layer2 = view.CreateLayer(outRaster,bands=[0,1,2])
end

