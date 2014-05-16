; docformat = 'rst'
; mjc_cursor_motion_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details

;+
; NAME:
;       MJC_CURSOR_MOTION
; PURPOSE:
;       Cursor communication with ENVI image windows
; AUTHOR;
;       Mort Canty (2006)
;       Juelich Research Center
;       m.canty@fz-juelich.de
; CALLING SEQUENCE:
;       MJC_Cursor_Motion, dn, xloc, yloc, $
;       xstart=xstart, ystart=ystart, event=event
; ARGUMENTS:
;       dn: display number
;       xloc,yloc: mouse position
; KEYWORDS
;       xstart, ystart: display origin
;       event: mouse event
; COMMON BLOCKS:
;       Cursor_Motion_C,dn,Cbtext,Rbtext,Cttext,Rttext
; DEPENDENCIES:
;       ENVI
;----------------------------------------------------
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

;+
; :Description:
;      Cursor communication with ENVI image windows
; :Params:
;      dn: in, required
;         display number
;      xloc,yloc: out, required
;         mouse position           
; :Keywords:
;      xstart, ystart: out, optional
;         display origin
;      event: in, required 
;         mouse event
; :Uses:
;      ENVI        
; :Author:
;      Mort Canty (2009) 
;-
pro mjc_cursor_motion, dn1, xloc, yloc, xstart=xstart, ystart=ystart, event=event
common cursor_motion_c, dn, Cbtext, Rbtext, Cttext, Rttext, Xtext, Ytext
  dn = dn1
; CalcHeight active?
  if (n_elements(Cbtext) eq 0) then Cbtext = -1L
  if (widget_info(Cbtext, /valid_id) eq 1) then begin
     if (event.press eq 1) and (event.modifiers eq 0) then begin
        widget_control, Cbtext, set_value = strtrim(xloc+1,2)
        widget_control, Rbtext, set_value = strtrim(yloc+1,2)
     endif
     if (event.press eq 1) and (event.modifiers eq 1) then begin
        widget_control, Cttext, set_value = strtrim(xloc+1,2)
        widget_control, Rttext, set_value = strtrim(yloc+1,2)
     endif
  endif
;MAD_scan_view, TCIMF_run or VIEWSEG_run active?
  if (n_elements(Xtext) eq 0) then Xtext = -1L
  if (widget_info(Xtext, /valid_id) eq 1) then begin
     envi_disp_query, dn, fid=fid
     if (event.press eq 1) and (event.modifiers eq 0) then begin
; place in text widgets
        widget_control, Xtext, set_value = strtrim(xstart+xloc+1,2)
        widget_control, Ytext, set_value = strtrim(ystart+yloc+1,2)
     endif
  end
end
