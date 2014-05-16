; docformat = 'rst'
; polsar__define.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.    
;+
; :Description:
;     Object class to store fully polarized SAR data in multilook covariance matrix form.
;     Property PLANES (sum of real and complex planes) identifies polarization type:
;        9: full polarimetry (upper diagonal part of 3x3 matrix: 3 real bands, 3 complex bands)
;        5: full polarimetry, azimuthal symmetry
;        3: full polarimetry, diagonal only
;        4: dual polarimetry
;        2: dual polarimetry, diagonal only
;        1: one-channel data      
; :Params:
;     fname: in, required
;        ENVI-style file header with samples, lines, bands, 
;           data type, interleave (BSQ) 
; :KEYWORDS:
;     NONE
; :Uses:
;     ENVI
; :Author:
;       Mort Canty (2013)        
;-
Function POLSAR::Init, fname, endian=endian
   compile_opt strictarr
; reads fully qualified ENVI header file name 
; for 6-band POLSAR covariance image and loads the data
   if n_elements(endian) eq 0 then endian = 'little' $
                              else endian = endian
   openR, lun, fname, /get_lun 
   aLine = ''
   readF, lun, aLine 
   if not (aLine eq 'ENVI') then begin
;    not an ENVI standard file, so assume DLR format
;     (mixed real and complex bands, root+polarization filename) 
      self.map_info = ptr_new(0)                                             
      i = strpos(aLine, '=')  
      self.samples = long(strmid(aline,i+1))
      readF, lun, aLine
      i = strpos(aLine, '=')  
      self.lines = long(strmid(aline,i+1))   
      readF, lun, aLine
      i = strpos(aLine, '=')  
      self.bands = long(strmid(aline,i+1))      
      free_lun, lun
      fname = file_dirname(fname)+'\'+file_basename(fname,'.hdr')   
      self.planes = 0  
      fn = fname+'C11'
      if file_test(fn) then begin
         data = read_binary(fn,data_type=4,data_dims=[self.samples,self.lines],endian=endian)
         if endian eq 'big' then data = reverse(transpose(data),2)
         self.planes++
      end else data = fltarr(self.samples,self.lines)
      self.C11 = ptr_new(data) 
      fn = fname+'C12'
      if file_test(fn) then begin
         data = read_binary(fn,data_type=6,data_dims=[self.samples,self.lines],endian=endian)
         if endian eq 'big' then data = reverse(transpose(data),2)         
         self.planes++
         self.planes++
      end else data = fltarr(self.samples,self.lines)
      self.C12 = ptr_new(data) 
      fn = fname+'hhvv'
      if file_test(fn) then begin
         data = read_binary(fn,data_type=6,data_dims=[self.samples,self.lines],endian=endian)
         if endian eq 'big' then data = reverse(transpose(data),2) 
         self.planes++
         self.planes++
      end else data = fltarr(self.samples,self.lines)
      self.C13 = ptr_new(data)  
      fn = fname+'hvhv'
      if file_test(fn) then begin
         data = read_binary(fn,data_type=4,data_dims=[self.samples,self.lines],endian=endian) 
         if endian eq 'big' then data = reverse(transpose(data),2)         
         self.planes++
      end else data = fltarr(self.samples,self.lines)
      self.C22 = ptr_new(data) 
      fn = fname+'hvvv'
      if file_test(fn) then begin
         data = read_binary(fn,data_type=6,data_dims=[self.samples,self.lines],endian=endian)
         if endian eq 'big' then data = reverse(transpose(data),2)        
         self.planes++
         self.planes++
      end else data = fltarr(self.samples,self.lines)
      self.C23 = ptr_new(data) 
      fn = fname+'vvvv'
      if file_test(fn) then begin
         data = read_binary(fn,data_type=4,data_dims=[self.samples,self.lines],endian=endian) 
         if endian eq 'big' then data = reverse(transpose(data),2)         
         self.planes++
      end else data = fltarr(self.samples,self.lines)
      self.C33 = ptr_new(data)
   end else begin
;    an ENVI standard file, complex data    
      free_lun, lun
      fname = file_dirname(fname)+'\'+file_basename(fname,'.hdr') 
      envi_open_file, fname, /no_interactive_query, /invisible, r_fid = fid
      envi_file_query, fid, nl=nl,ns=ns,nb=nb,dims=dims,data_type=dt
      self.map_info = ptr_new(envi_get_map_info(fid=fid))
      self.lines=nl
      self.samples=ns
      self.bands=nb
      planes = 0
      if (nb ne 6) or ((dt ne 6) and (dt ne 9)) then  return, 0
      data = real_part(envi_get_data(/complex, fid=fid, dims=dims, pos=0))
      self.C11 = ptr_new(data)
      if total(data) ne 0 then planes++
      data = envi_get_data(/complex, fid=fid, dims=dims, pos=1)
      self.C12 = ptr_new(data)
      if total(data) ne 0 then planes=planes+2
      data = envi_get_data(/complex, fid=fid, dims=dims, pos=2)
      self.C13 = ptr_new(data)
      if total(data) ne 0 then planes=planes+2
      data = real_part(envi_get_data(/complex, fid=fid, dims=dims, pos=3))
      self.C22 = ptr_new(data)
      if total(data) ne 0 then planes++
      data = envi_get_data(/complex, fid=fid, dims=dims, pos=4)
      self.C23 = ptr_new(data)
      if total(data) ne 0 then planes=planes+2 
      data = real_part(envi_get_data(/complex, fid=fid, dims=dims, pos=5))
      self.C33 = ptr_new(data)
      if total(data) ne 0 then planes++ 
      self.planes = planes                 
   endelse 
   case self.planes of
     9: self.d = 3
     5: self.d = 3
     3: self.d = 3
     4: self.d = 2
     2: self.d = 2
     1: self.d = 1
   endcase    
   return, 1 
End

Pro POLSAR::getProperty, lines=lines, samples=samples, bands=bands,d=d, $
    planes=planes, span=span, c11=c11, c22=c22, c33=c33, map_info=map_info 
   compile_opt strictarr
   
   if arg_present(lines) then lines=self.lines
   if arg_present(samples) then samples=self.samples
   if arg_present(bands) then bands=self.bands
   if arg_present(d) then d=self.d
   if arg_present(planes) then planes=self.planes
   if arg_present(span) then span = *(self.c11)+2*(*self.c22)+*(self.c33)
   if arg_present(c11) then c11 = complex(*(self.c11))
   if arg_present(c22) then c22 = complex(*(self.c22))
   if arg_present(c33) then c33 = complex(*(self.c33))
   if arg_present(map_info) then map_info = *(self.map_info)
End   

Pro POLSAR::ToENVI, fid = r_fid
   compile_opt strictarr

   bnames = ['C11','C12','C13','C22','C23','C33']
   im = complexarr(self.samples,self.lines,6)
   im[*,*,0] = *(self.c11)
   im[*,*,1] = *(self.c12)
   im[*,*,2] = *(self.c13)
   im[*,*,3] = *(self.c22)
   im[*,*,4] = *(self.c23)
   im[*,*,5] = *(self.c33) 
   tmp = reform(im,self.lines*self.samples,6)
   idx = where(total(tmp,1) eq 0,count)
   if count gt 0 then bnames[idx] = '-'  
   ENVI_Enter_Data, im, bnames=bnames, r_fid=r_fid    
End

Pro POLSAR::CleanUp
   ptr_free, self.c11
   ptr_free, self.c12
   ptr_free, self.c13
   ptr_free, self.c22
   ptr_free, self.c23
   ptr_free, self.c33
End   

Pro POLSAR::times, n
; multiply covariance matrix by n
;  (if n is ENL, result is complex Wishart distributed)
   compile_opt strictarr
   *(self.C11) = *(self.C11)*n
   *(self.C12) = *(self.C12)*n
   *(self.C13) = *(self.C13)*n
   *(self.C22) = *(self.C22)*n
   *(self.C23) = *(self.C23)*n
   *(self.C33) = *(self.C33)*n
End   

Function POLSAR::Band, band
   compile_opt strictarr
   case band of
      'C11': return, *(self.C11)
      'C12': return, *(self.C12)
      'C13': return, *(self.C13)
      'C22': return, *(self.C22)
      'C23': return, *(self.C23) 
      'C33': return, *(self.C33)
   endcase
End         

Pro POLSAR__Define
compile_opt strictarr
class =  { POLSAR,          $
           samples: 0L,     $
           lines:   0L,     $
           bands:   0L,     $
           planes:  0L,     $
           d:       0L,     $ 
           map_info: Ptr_New(), $
           C11: Ptr_New(), $
           C12: Ptr_New(), $
           C13: Ptr_New(), $
           C22: Ptr_New(), $
           C23: Ptr_New(), $
           C33: Ptr_New() }
End