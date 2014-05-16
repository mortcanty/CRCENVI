;+
; NAME:
;       WARP_SHIFT
; PURPOSE:
;       Use RST with bilinear interpolation to shift band to sub-pixel accuracy
; AUTHOR
;       Mort Canty (2004)
;       Juelich Research Center
;       m.canty@fz-juelich.de
; CALLING SEQUENCE:
;       sband = Warp_Shift(band,shft)
; ARGUMENTS:
;       band: the image band to be shifted
;       shft: array of length 2 giving the amount to be shifted
; KEYWORDS:
;       None
; DEPENDENCIES:
;       ENVI
;---------------------------------------------------------------------------
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program; if not, write to the Free Software
;    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
;----------------------------------------------------------

function warp_shift, band, shft
; if integer shift, use SHIFT function and no resampling
   if (fix(shft[0]) eq shft[0]) and (fix(shft[1]) eq shft[1]) then return, shift(band,shft)
; otherwise use ENVI batch program
   sz = size(band)
   dims = [-1L,0,sz[1]-1,0,sz[2]-1]
   pts = dblarr(4,4)
   pts[0,0] =   sz[1]/4   & pts[1,0] =   sz[2]/4
   pts[0,1] = 3*sz[1]/4   & pts[1,1] =   sz[2]/4
   pts[0,2] =   sz[1]/4   & pts[1,2] = 3*sz[2]/4
   pts[0,3] = 3*sz[1]/4   & pts[1,3] = 3*sz[2]/4
   for i=0,3 do pts[2:3,i] = pts[0:1,i] - shft
   envi_enter_data, band, r_fid=r_fid
   envi_doit, 'envi_register_doit', b_fid=r_fid, w_fid=r_fid, $
               w_dims=dims, w_pos=0, /in_memory, pts=pts,  $
               method=2, $ ; RST with cubic convolution
               r_fid=r_fid1, x0=0, y0=0, $
               xsize = sz[1], ysize = sz[2]
   envi_file_mng, id=r_fid, /remove
   result = envi_get_data(fid=r_fid1, dims=dims, pos=0)
   envi_file_mng, id=r_fid1, /remove
   return, result
end


