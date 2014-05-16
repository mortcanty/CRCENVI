; docformat = 'rst'
; ENLsim.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
    
pro ENLsim
; utility to simulate a polSAR image in covariance matrix format
; Mort Canty, with help from Knut Conradsen (2013)

   COMPILE_OPT IDL2
   
;  envi_select, title='Choose (spatial subset of) representative quad polSAR image', $
;                               fid=fid, pos=pos, dims=dims
;                                              
;  if (fid eq -1) then begin
;      print, 'cancelled'
;     return
;  endif
;
;  envi_file_query, fid,fname=fname
;  cols = dims[2]-dims[1]+1
;  rows = dims[4]-dims[3]+1
;  bands = 6
;  pixels = cols*rows
;; data matrix
;  G = complexarr(bands,pixels)
;; read in image
;  for i=0,bands-1 do G[i,*] = envi_get_data(fid=fid,dims=dims,pos=pos[i],/complex)
;; determine average covariance matrix  
;  Sigma = complexarr(3,3)
;  Sigma[0,0] = mean(G[0,*])
;  Sigma[1,0] = mean(G[1,*])
;  Sigma[2,0] = mean(G[2,*])
;  Sigma[1,1] = mean(G[3,*])
;  Sigma[2,1] = mean(G[4,*])
;  Sigma[2,2] = mean(G[5,*])
;; replace Sigma by its lower Cholesky factorization   
;  LA_CHOLDC, Sigma,/upper
;  C = Sigma 
   
  n = 12
  p = 3
  seed = 12321
   
  outimage = complexarr(500^2,6)
   
  for i = 0L, 500^2-1 do begin

    X = randomu(seed,p,n,/normal)
    Y = randomu(seed,p,n,/normal)
    
    S = complex( transpose(X)##X + transpose(Y)##Y, -(transpose(X)##Y - transpose(Y)##X) )/2.0
    
;    W = complex(transpose(C))##S##C 
    W = S
    
    case p of 
       1: outimage[i,0] = W[0]
       2: outimage[i,[0,1,3]] = [W[0,0],W[1,0],W[1,1]]
       3: outimage[i,*] = [W[0,0],W[1,0],W[2,0],W[1,1],W[2,1],W[2,2]] / n
    endcase
    
  endfor 
  
  case p of
     1: bnames = ['HHHH','-','-','-','-','-']
     2: bnames = ['HHHH','HHHV','-','HVHV','-','-']
     3: bnames = ['HHHH','HHHV','HHVV','HVHV','HVVV','VVVV']
  endcase   
  
  envi_enter_data, reform(outimage,500,500,6), bnames=bnames
   
end