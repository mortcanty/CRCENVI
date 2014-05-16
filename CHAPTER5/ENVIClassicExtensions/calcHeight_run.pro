; docformat = 'rst'
; calcheight_run.pro
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.

PRO calcHeight_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      /DISPLAY, $
      VALUE = 'Height Calculator', $
      REF_VALUE = 'Tools', $
      EVENT_PRO = 'calcHeight_run', $
      UVALUE = 'CALCHEIGHT',$
      POSITION = 'last', $
      /SEPARATOR
END

Function parseRPBline, aLine
   pStart = strpos(aLine,'=')
   pEnd = strpos(aLine,';')
   if pEnd eq -1 then pEnd = strpos(aLine,',')
   r = 0.0D
   readS, strmid(aLine,pStart+1,pEnd-pStart-1), r
   return, r
end

Function parseRPCline, aLine
   pStart = strpos(aLine,': ')
   pEnd = strpos(aLine,'pixels')
   if pEnd eq -1 then pEnd = strpos(aLine,'degrees')
   if pEnd eq -1 then pEnd = strpos(aLine,'meters')
   if pEnd eq -1 then pEnd = strlen(aLine)
   r = 0.0D
   readS, strmid(aLine,pStart+1,pEnd-pStart-1), r
   return, r
end

function parseRPB, fileName
   !error_state.code=0
   catch, theError
   if theError ne 0 then return, {errBias: -100}
   openR, lun, fileName, /get_lun
; skip the first 4 lines
   header = strarr(4)
   readF, lun, header
   aLine = ''
   readF, lun, aLine
   errBias = parseRPBline(aLine)
   readF, lun, aLine
   errRand = parseRPBline(aLine)
   readF, lun, aLine
   lineOffset = parseRPBline(aLine)
   readF, lun, aLine
   sampOffset = parseRPBline(aLine)
   readF, lun, aLine
   latOffset = parseRPBline(aLine)
   readF, lun, aLine
   longOffset = parseRPBline(aLine)
   readF, lun, aLine
   heightOffset = parseRPBline(aLine)
   readF, lun, aLine
   lineScale = parseRPBline(aLine)
   readF, lun, aLine
   sampScale = parseRPBline(aLine)
   readF, lun, aLine
   latScale = parseRPBline(aLine)
   readF, lun, aLine
   longScale = parseRPBline(aLine)
   readF, lun, aLine
   heightScale = parseRPBline(aLine)
   readF, lun, aLine
   lineNumCoef =dblarr(20)
   for i=0,19 do begin
      readF, lun, aLine
      lineNumCoef[i] = parseRPBline(aLine)
   endfor
   readF,lun,aLine
   lineDenCoef =dblarr(20)
   for i=0,19 do begin
      readF, lun, aLine
      lineDenCoef[i] = parseRPBline(aLine)
   endfor
   readF,lun,aLine
   sampNumCoef =dblarr(20)
   for i=0,19 do begin
      readF, lun, aLine
      sampNumCoef[i] = parseRPBline(aLine)
   endfor
   readF,lun,aLine
   sampDenCoef =dblarr(20)
   for i=0,19 do begin
      readF, lun, aLine
      sampDenCoef[i] = parseRPBline(aLine)
   endfor
   free_lun, lun
   return,   {errBias:errBias,errRand:errRand,lineOffset:lineOffset,sampOffset:sampOffset, $
              latOffset:latOffset,longOffset:longOffset,heightOffset:heightOffset,$
              lineScale:lineScale,sampScale:sampScale, $
              latScale:latScale,longScale:longScale,heightScale:heightScale,$
              lineNumCoef:lineNumCoef,lineDenCoef:lineDenCoef,$
              sampNumCoef:sampNumCoef,sampDenCoef:sampDenCoef}
end

function parseRPC, fileName
   !error_state.code=0
   catch, theError
   if theError ne 0 then return, {errBias:-100}
   openR, lun, fileName, /get_lun
; skip the first 4 lines
   aLine = ''
   readF, lun, aLine
   lineOffset = parseRPCline(aLine)
   readF, lun, aLine
   sampOffset = parseRPCline(aLine)
   readF, lun, aLine
   latOffset = parseRPCline(aLine)
   readF, lun, aLine
   longOffset = parseRPCline(aLine)
   readF, lun, aLine
   heightOffset = parseRPCline(aLine)
   readF, lun, aLine
   lineScale = parseRPCline(aLine)
   readF, lun, aLine
   sampScale = parseRPCline(aLine)
   readF, lun, aLine
   latScale = parseRPCline(aLine)
   readF, lun, aLine
   longScale = parseRPCline(aLine)
   readF, lun, aLine
   heightScale = parseRPCline(aLine)
   lineNumCoef =dblarr(20)
   for i=0,19 do begin
      readF, lun, aLine
      lineNumCoef[i] = parseRPCline(aLine)
   endfor
   lineDenCoef =dblarr(20)
   for i=0,19 do begin
      readF, lun, aLine
      lineDenCoef[i] = parseRPCline(aLine)
   endfor
   sampNumCoef =dblarr(20)
   for i=0,19 do begin
      readF, lun, aLine
      sampNumCoef[i] = parseRPCline(aLine)
   endfor
   sampDenCoef =dblarr(20)
   for i=0,19 do begin
      readF, lun, aLine
      sampDenCoef[i] = parseRPCline(aLine)
   endfor
   free_lun, lun
   return,   {errBias:0,errRand:0,lineOffset:lineOffset,sampOffset:sampOffset, $
              latOffset:latOffset,longOffset:longOffset,heightOffset:heightOffset,$
              lineScale:lineScale,sampScale:sampScale, $
              latScale:latScale,longScale:longScale,heightScale:heightScale,$
              lineNumCoef:lineNumCoef,lineDenCoef:lineDenCoef,$
              sampNumCoef:sampNumCoef,sampDenCoef:sampDenCoef}
end

; routines for Newton's method
;--------------------------------------------------------------

function polynom, V, a
; return a cubic polynomial in X (long) Y (lat) Z (height) with coefficients a[]
; See Grodecki, Ikonos stereo feature extraction - RPC approach, p2
   X=V[0]
   Y=V[1]
   Z=V[2]
   return, a[0]+a[1]*X+a[2]*Y+a[3]*Z+a[4]*X*Y+a[5]*X*Z+a[6]*Y*Z+a[7]*X*X+a[8]*Y*Y+a[9]*Z*Z+$
         a[10]*X*Y*Z+a[11]*X*X*X+a[12]*X*Y*Y+a[13]*X*Z*Z+a[14]*X*X*Y+a[15]*Y*Y*Y+a[16]*Y*Z*Z+$
         a[17]*X*X*Z+a[18]*Y*Y*Z+a[19]*Z*Z*Z
end

function newtFunc, X
common shared, RFM, Cb, Rb, Ct, Rt, elev
; unknowns are X[0]=Long, X[1]=Lat, Height is fixed at elev
; Cb and Rb are scaled
   a = RFM.lineNumCoef
   b = RFM.lineDenCoef
   c = RFM.sampNumCoef
   d = RFM.sampDenCoef
   h = (elev-RFM.heightOffset)/RFM.heightScale
; evaluate RPs for row (f) and column (g)
   f = polynom([X[0],X[1],H],a)/polynom([X[0],X[1],H],b)
   g = polynom([X[0],X[1],H],c)/polynom([X[0],X[1],H],d)
   return, [f-Rb,g-Cb]
end

function height, X
common shared, RFM, Cb, Rb, Ct, Rt, elev
; X[0]=Long, X[1]=Lat, X[2]=height
; Ct and Rt are absolute pixel values
   a = RFM.lineNumCoef
   b = RFM.lineDenCoef
   c = RFM.sampNumCoef
   d = RFM.sampDenCoef
; evaluate RPs for row (f) and column (g)
   f = polynom(X,a)/polynom(X,b)
   g = polynom(X,c)/polynom(X,d)
   R = RFM.lineScale*f+RFM.lineOffset
   C = RFM.sampScale*g+RFM.sampOffset
   return, sqrt((R-Rt)^2+(C-Ct)^2)
end

; Event handlers for GUI
;--------------------------------------------------------

PRO calcHeight_cleanup, tlb
   widget_control,tlb,Get_Uvalue=state,/no_copy
   if n_elements(state) eq 0 then return
   ptr_free, state.RFM
END

PRO onQuit, event
   widget_control,event.top,Get_Uvalue=state,/no_copy
   if n_elements(state) eq 0 then return
   ptr_free, state.RFM
   widget_control,event.top,/destroy
END

PRO onLoadRPC, event
   widget_control,event.top,Get_Uvalue=state,/no_copy
   fileName = Dialog_Pickfile(Filter=['*.rpb','*.rpc'],/Read)
   if fileName eq '' then begin
      widget_control,event.top,Set_Uvalue=state,/no_copy
      return
   endif
   suffix = ''
   readS, strmid(fileName,strlen(fileName)-3,strlen(fileName)-1), suffix
   if (suffix eq 'rpc') or (suffix eq 'RPC') then result = parseRPC(fileName) $
                                             else result = parseRPB(fileName)
   if result.errBias eq -100 then begin
       void = dialog_message('Error reading RPC file',/error)
       widget_control,event.top,Set_Uvalue=state,/no_copy
       return
   endif
   if result.lineNumCoef[3] eq 0 then begin
       void = dialog_message('Image is orthorectified',/error)
       widget_control,event.top,Set_Uvalue=state,/no_copy
       return
   endif

   state.RFM = ptr_new(result)
   widget_control,state.calcID,sensitive=1
   widget_control,state.textID,set_value=fileName
   widget_control,event.top,Set_Uvalue=state,/no_copy
END

PRO onLoadDEM, event
   widget_control,event.top,Get_Uvalue=state,/no_copy
   envi_select, title='Choose DEM file', fid=DEMfid, /band
   if (DEMfid eq -1) then begin
      widget_control,event.top,Set_Uvalue=state,/no_copy
      return
   endif
   state.DEMfid=DEMfid
   widget_control,event.top,Set_Uvalue=state,/no_copy
END

PRO onHelp, event
   XDISPLAYFILE, 'calcheight_help.txt',group=event.top, height=20, width=50, $
   text=[ $
'Help on CalcHeight:', $
'===================', $
' ', $
'After loading an RPC file, click on the bottom of', $
'a vertical structure to set the base height and then', $
'shift-click on the top of the structure. Press', $
'the CALC button to display the structures height', $
'latitude, longitude and base elevation. The number', $
'in brackets next to the height is the minimum distance', $
'(in pixels) between the top pixel and a vertical line', $
'through the bottom pixel. It should be of the order of 1', $
'or less.', $
' ', $
'If no DEM is loaded, the base elevation is the average', $
'value for the whole scene.', $
' ', $
'If a DEM is used, the base elevation is taken from it.', $
'The latitude and longitude are then orthorectified values.']
END

PRO onAbout, event
   dummy=dialog_message(['CalcHeight ',$
                         'M. Canty, 2005'],/info)
END

PRO onCalc, event
common shared, RFM, Cb, Rb, Ct, Rt, elev
common cursor_motion_c, dn, Cbtext, Rbtext, Cttext, Rttext, Xtext, Ytext
   !error_state.code=0
   catch,error
   if error NE 0 then begin
     help,/last_message,output=traceback
     errarray=['Error Caught',traceback]
     void = dialog_message(errarray,/error)
     return
   endif
   widget_control,event.top,get_uvalue=state,/no_copy
; obtain pixel values from text widgets:
   widget_control,state.Cbtext,get_value=Cb
   widget_control,state.Rbtext,get_value=Rb
   widget_control,state.Cttext,get_value=Ct
   widget_control,state.Rttext,get_value=Rt
; convert from string array to double-precision scalar:
   Cb=double(Cb[0])
   Rb=double(Rb[0])
   Ct=double(Ct[0])
   Rt=double(Rt[0])
; set RFM in common block
   RFM = *state.RFM
; set local elevation in common block
   if state.DEMfid eq -1 then elev = RFM.heightOffset $
   else begin
   ; get the map coordinates of the base
      envi_disp_query, dn, fid=fid
      envi_convert_file_coordinates, fid[0], Cb, Rb, e, n, /to_map
   ; get corresonding pixels coordinates in the DEM
      envi_convert_file_coordinates, state.DEMfid, X, Y, e, n
   ; extract the absolute elevation from the DEM
      X=round(X)
      Y=round(Y)
      dims = [-1,X,X,Y,Y]
      elev = (envi_get_data(fid=state.DEMfid,dims=dims,pos=0))[0]
   endelse
; set scaled pixel values for Cb and Rb in common block
   Cb = (Cb-RFM.sampOffset)/RFM.sampScale
   Rb = (Rb-RFM.lineOffset)/RFM.lineScale
; estimates for [Xb,Yb]
   X = dblarr(2)
   X = NEWTON(X,'newtFunc',check=check,/double)
; get best matching height
   matches = fltarr(300)
   for h=0,299 do begin
      Z = (h + elev - RFM.heightOffset)/RFM.heightScale
      matches[h] = height([X[0],X[1],Z])
   endfor
   h_best = (where(matches eq min(matches)))[0]
 ; set result
   latitude = X[1]*(*state.RFM).latScale + (*state.RFM).latOffset
   longitude = X[0]*(*state.RFM).longScale + (*state.RFM).longOffset
   s = string(matches[h_best],format='(F5.2)')
   widget_control,state.heighttext,set_value=STRTRIM(h_best,2)+' m (' + s + ' pxl)'
   widget_control,state.lattext,set_value=STRTRIM(latitude,2)
   widget_control,state.longtext,set_value=STRTRIM(longitude,2)
   widget_control,state.elevtext,set_value=STRTRIM(round(elev),2)+' m'
   widget_control,event.top,Set_Uvalue=state,/no_copy
END

PRO NULL_EVENTS,event
END

;+
; :Description:
;      ENVI extention to determine height of vertical
;      structures in high resolution images using RFMs
; :Params:
;      event:  in, required 
;         if called from the ENVI menu  
; :Uses:
;      ENVI, MJC_CURSOR_MOTION
; :Author:
;      Mort Canty (2013)         
;-
PRO calcHeight_run,event

; cursor communication
common cursor_motion_c, dn, Cbtext, Rbtext, Cttext, Rttext

  tlb=widget_base(title='CalcHeight',/column,mbar=menubarID,xoffset=100,yoffset=100)
  fileID        = widget_button(menubarID,value='File',/menu)
  helpTopID     = widget_button(menubarID,value='Help',/menu)
  helpID        = widget_button (helpTopID, value='Help', event_pro='ONHelp')
  aboutID       = widget_button (helpTopID, value='About', event_pro='ONAbout')
  loadRPCID     = widget_button (fileID, value='Load RPC file', event_pro='ONloadRPC')
  loadDEMID     = widget_button (fileID, value='Load DEM file', event_pro='ONloadDEM')
  quitID        = widget_button(fileID,value='Quit',event_pro='ONquit')
  textID        = widget_text(tlb,value='No RPC file loaded')
  row1=widget_base(tlb,/row)
    label1=widget_label(row1,VALUE='X-value of base: ')
    Cbtext=widget_text(row1,xsize=10,ysize=1,EVENT_PRO='NULL_EVENTS',/editable)
    label2=widget_label(row1,VALUE=' Y-value of base: ')
    Rbtext=widget_text(row1,xsize=10,ysize=1,EVENT_PRO='NULL_EVENTS',/editable)
  row2=widget_base(tlb,/row)
    label3=widget_label(row2,VALUE='X-value of top:     ')
    Cttext=widget_text(row2,xsize=10,ysize=1,EVENT_PRO='NULL_EVENTS',/editable)
    label4=widget_label(row2,VALUE=' Y-value of top:     ')
    Rttext=widget_text(row2,xsize=10,ysize=1,EVENT_PRO='NULL_EVENTS',/editable)
  row3=widget_base(tlb,/row)
    label5=widget_label(row3,VALUE='Height: ')
    heighttext=widget_text(row3,xsize=12,ysize=1,EVENT_PRO='NULL_EVENTS')
    label6=widget_label(row3,VALUE=' Lat: ')
    lattext=widget_text(row3,xsize=10,ysize=1,EVENT_PRO='NULL_EVENTS')
    label6=widget_label(row3,VALUE=' Long: ')
    longtext=widget_text(row3,xsize=10,ysize=1,EVENT_PRO='NULL_EVENTS')
    label6=widget_label(row3,VALUE=' Elev: ')
    elevtext=widget_text(row3,xsize=10,ysize=1,EVENT_PRO='NULL_EVENTS')
  row6=widget_base(tlb,/row)
    calcID=widget_button(row6,VALUE='Compute',EVENT_PRO='ONCalc',sensitive=0)

    RFM = {dummy:1}
    DEMfid = -1

; realize the GUI
   widget_control,tlb,/realize

; initialize state structure:
   state={Cbtext:Cbtext,Rbtext:Rbtext,Cttext:Cttext,Rttext:Rttext, $
          HeightText:HeightText,LatText:LatText,LongText:LongText,ElevText:ElevText, $
          RFM:ptr_new(),calcID:calcID,textID:textID,DEMfid:DEMfid}

; use state structure as widget uvalues:
   widget_control,tlb,set_uvalue=state,group_leader=event.top

; start XMANAGER:
   XMANAGER,'CalcHeight',tlb, /no_block, cleanup='CalcHeight_cleanup'

END



