pro ex1_2
COMPILE_OPT IDL2

   fname = dialog_pickfile(filter='*.xml', $
            title='Select COSAR image')
   if fname eq '' then return
   tdx_openfile, /no_envi, filename=fname, $
                              saInfo=saInfo
   data = tdx_readspatial(band =[0],$
                          xs=1000, $
                          xe=5999, $
                          ys=1000,$
                          ye=5999,$ 
                          saInfo=saInfo)					  
   tdx_closefile, saInfo    
   window, 11, xsize=5000, ysize=5000  
   tvscl, abs(data)                            

end