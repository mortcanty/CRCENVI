;docformat = 'rst'
; aaa_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO aaa_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'GPULib API', $
      REF_VALUE = 'About ENVI', $
      EVENT_PRO = 'aaa_run', $
      UVALUE = 'AAA',$
      POSITION = 'after' 
END

;+
; :Description:
;    Place in SAVE_ADD to set
;    a menu item under ENVI Help which opens
;    a link to the GPULib API using
;    Mike Galloy's mg_open_url procedure
; :Params:
;       event:  in, optional 
;          required if called from ENVI                 
; :Uses:
;       ENVI::
;       MG_OPEN_URL::       
; :Author:
;       Mort Canty (2009)      
;-
pro aaa_run,event
   mg_open_url,'http://www.txcorp.com/products/GPULib/idl_docs/index.html'
end 