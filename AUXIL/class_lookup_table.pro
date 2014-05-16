; docformat = 'rst'
; class_lookup_table.pro
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
;       returns colors from a 16-color lookup table     
; :Params:
;       ptr:  in, required 
;          color indices          
; :Author:
;       Mort Canty (2009)    
;-
Function class_lookup_table, ptr
   black =     [0,0,0]
   white =     [255,255,255]
   red   =     [255,0,0]
   green =     [0,255,0]
   blue  =     [0,0,255]
   yellow=     [255,255,0]
   cyan  =     [0,255,255]
   magenta =   [255,0,255]
   maroon =    [176,48,96]
   seagreen =  [46,139,87]
   purple =    [160,32,240]
   coral =     [255,127,80]
   aquamarine =[127,255,212]
   orchid =    [218,112,214]
   sienna =    [160,82,45]
   chartreuse =[127,255,0]
   table = [[black],[white],[red],[green],[blue],[yellow],[cyan],[magenta],$
             [maroon],[seagreen],[purple],[coral],[aquamarine],[orchid],$
             [sienna],[chartreuse]]
   return, table[*,ptr]
End