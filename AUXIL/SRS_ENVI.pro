;  Authors:
;  Hongjie Xie,  Nigel Hicks,  G. Randy Keller, Haitao Huang, Vladik Kreinovich
;
;  Title of the paper:
;  An IDL/ENVI implementation of the FFT Based Algorithm for Automatic Image Registration
;
;
; PURPOSE:
;       The purpose of this program is to provide a quick-and-easy way
;       to automatic register images.
;
; CATEGORY:
;       ENVI User Function, Visual / Display, Map Information
;
; CALLING SEQUENCE:
;       SRS_ENVI
;
; INPUT PARAMETERS:
;       ENVI event structure.
;
;
; OUTPUTS:
;       ENVI image file or array in memory.
;
;
;
; MODIFICATION HISTORY:
;       Written IDL/ENVI Team at UTEP in July 2000, updated in May 2002.
;-

PRO SRS_ENVI_event,event
;main widget event handler:
;error checking:
!error_state.code=0
catch,error
if error NE 0 then begin
  help,/last_message,output=traceback
  errarray=['Error Caught',traceback]
  dummy=dialog_message(errarray,/error,/cancel)
  if STRUPCASE(dummy) EQ 'CANCEL' then return
endif
;get state structure:
widget_control,event.top,get_uvalue=state
;check to see if TLB was killed:
if TAG_NAMES(event,/STRUCTURE_NAME) EQ 'WIDGET_KILL_REQUEST' then begin
  device,decomposed=state.olddc
  widget_control,event.top,/destroy
  return
endif
;event handler for widget_outfm and colorbar widget_draw expose events:
case (event.id) of
  state.fileio:    BEGIN ;widget_outfm event:
                     if event.result.in_memory EQ 0 then begin ;output to file:
                       state.saveflag=0B
                       widget_control,state.fileio,get_value=current
                       if (event.result.name EQ '') AND (current.name EQ '') then begin
                         dummy=dialog_message(['You have to enter an output',$
                                               'filename before running SRS_ENVI !'])
                         widget_control,event.top,set_uvalue=state
                         return
                       endif
                       outfind=FINDFILE(event.result.name)
                       outexist=size(outfind,/n_dimensions)
                       if outexist EQ 1 then begin
                         dummy=dialog_message([event.result.name,'This file already exists.','',$
                                              'Replace existing file?'],/question,/default_no)
                         if strupcase(dummy) EQ 'NO' then begin
                           widget_control,state.fileio,set_value={out_name:'',in_memory:0L}
                           state.outfile=''
                           widget_control,event.top,set_uvalue=state
                           return
                         endif else if strupcase(dummy) EQ 'YES' then begin
                           state.outfile=event.result.name
                         endif
                       endif else begin
                         state.outfile=event.result.name
                       endelse
                     endif else begin ;output to memory
                       state.saveflag=1B
                     endelse
                   END

endcase
widget_control,event.top,set_uvalue=state
END


PRO SRS_ENVIOpenBaseim,event
;event handler for Bim "Import Band..." button:
widget_control,event.top,get_uvalue=state
envi_select,title='Select Base Image',fid=fid,pos=pos,dims=dims,/no_dims,/band_only
envi_file_query,fid,fname=fname,ns=ns,nl=nl,pixel_size=ps,h_map=h_map
widget_control,state.text1,set_value=fname
map={envi_map_struct}
if h_map NE 0 then begin
  handle_value,h_map,map
  state.mapflag=1B
endif
state.Bimdims=dims
state.Bimfid=fid
state.Bimname=fname
state.Bimpos=pos
state.Bimns=ns
state.Bimnl=nl
state.Bimmap=map
if ps[0] NE 0 and ps[1] NE 0 then state.psflag=1B
if ps[0] EQ 0 and h_map NE 0 then begin
  state.Bimps=map.ps
endif else begin
  state.Bimps=ps
endelse
if state.previewflag EQ 1B then begin
  state.previewflag=0B
  state.preview_id=0L
  widget_control,state.tlb,set_uvalue=state
  widget_control,state.tlb5,/destroy
endif
state.stretchflag=0B
widget_control,event.top,set_uvalue=state
END

PRO SRS_ENVIRemoveBaseim,event
;event handler for Bim "Remove" button:
widget_control,event.top,get_uvalue=state
if state.previewflag EQ 1B then begin
  state.previewflag=0B
  state.preview_id=0L
  widget_control,state.tlb,set_uvalue=state
  widget_control,state.tlb5,/destroy
endif
widget_control,state.text1,set_value=''
state.Bimfid=-1L
state.stretchflag=0B
state.mapflag=0B
state.psflag=0B
widget_control,event.top,set_uvalue=state
END


PRO OpenSRS_ENVIim,event
;event handler for image "Import Band..." button:
widget_control,event.top,get_uvalue=state
envi_select,title='Select Image with scaling, rotation and shift',fid=fid,pos=pos,dims=dims,/no_dims,/band_only
envi_file_query,fid,fname=fname,ns=ns,nl=nl,pixel_size=ps,h_map=h_map
widget_control,state.text2,set_value=fname
map={envi_map_struct}
if h_map NE 0 then begin
  handle_value,h_map,map
  state.mapflag=1B
endif
state.imagedims=dims
state.imagefid=fid
state.imagename=fname
state.imagepos=pos
state.imagens=ns
state.imagenl=nl
state.imagemap=map
if ps[0] NE 0 and ps[1] NE 0 then state.psflag=1B
if ps[0] EQ 0 and h_map NE 0 then begin
  state.imageps=map.ps
endif else begin
  state.imageps=ps
endelse
if state.previewflag EQ 1B then begin
  state.previewflag=0B
  state.preview_id=0L
  widget_control,state.tlb,set_uvalue=state
  widget_control,state.tlb5,/destroy
endif
state.stretchflag=0B
widget_control,event.top,set_uvalue=state
END


PRO RemoveSRS_ENVIim,event
;event handler for Bim "Remove" button:
widget_control,event.top,get_uvalue=state
if state.previewflag EQ 1B then begin
  state.previewflag=0B
  state.preview_id=0L
  widget_control,state.tlb,set_uvalue=state
  widget_control,state.tlb5,/destroy
endif
widget_control,state.text2,set_value=''
state.imagefid=-1L
state.stretchflag=0B
state.mapflag=0B
state.psflag=0B
widget_control,event.top,set_uvalue=state
end


FUNCTION Ritio, ARR1, ARR2
  ArrSize = SIZE(ARR1)
  Cols = ArrSize[1]
  Rows = ArrSize[2]

  R=COMPLEXARR(Cols,Rows,/Nozero)
  R1=COMPLEXARR(Cols,Rows,/Nozero)
  R2=FLTARR(Cols,Rows,/Nozero)

  R1=ARR1*(temporary(conj(ARR2)))
  R2=(temporary(abs(ARR1)))*(temporary(abs(ARR2)))

  for i = 0, Rows - 1 do begin
    for j = 0, Cols - 1 do begin
      JUNK = CHECK_MATH()
      !EXCEPT=0
      R[j, i] = temporary(R1[j, i]) / temporary(R2[j, i])
      if FINITE(temporary(FLOAT(R[j, i]))) EQ 0 then begin
      endif
      if FINITE(temporary(imaginary(R[j, i]))) EQ 0 then begin
      endif
    endfor
  endfor
  RETURN, R
END

FUNCTION LogPolar, RECT, LP
  Pi = 3.14159365359
  RArrSize = SIZE(RECT)
  RCols = RArrSize[1]
  RRows = RArrSize[2]
  LP=temporary(FLTARR(RCols, RRows,/Nozero))

dTheta= 1.0 * pi / RRows          ;the step of angle
b = 10 ^ (alog10(RCols) / RCols)  ;The base for the log-polar conversion from rectangle.
for i=0.0, RRows-1.0 do begin
  Theta=i * dTheta
  for j=0.0, RCols-1.0 do begin
    r = b ^ j - 1                 ;The log-polar
    x=r * cos(Theta) + RCols / 2.0
    y=r * sin(Theta) + RRows / 2.0
    x0=floor(x)
    y0=floor(y)
    x1=x0+1
    y1=y0+1
    if (x0 LE RCols-1) and (y0 LE RRows-1) and (x0 GT 1) and (y0 GT 1) THEN V00=RECT[x0,y0]$
                                                                       else V00=0.0
    if (x0 LE RCols-1) and (y1 LE RRows-1) and (x0 GT 1) and (y1 GT 1) THEN V01=RECT[x0,y1]$
                                                                       else V01=0.0
    if (x1 LE RCols-1) and (y0 LE RRows-1) and (x1 GT 1) and (y0 GT 1) THEN V10=RECT[x1,y0]$
                                                                       else V10=0.0
    if (x1 LE RCols-1) and (y1 LE RRows-1) and (x1 GT 1) and (y1 GT 1) THEN V11=RECT[x1,y1]$
                                                                       else V11=0.0
    V=V00*(x1-x)*(y1-y)+V01*(x1-x)*(y-y0)+V10*(x-x0)*(y1-y)+V11*(y-y0)*(x-x0);Bilinear interpolation
    LP[j,i]=temporary(V)
  endfor
endfor
RETURN, LP
END

PRO HighPass, ARR
  Pi = 3.14159365359
  ArrSize = SIZE(ARR)
  Cols = ArrSize[1]
  Rows = ArrSize[2]
  hpMatrix = temporary(FLTARR((Cols+1)*(Rows+1)))
  for i=0, Rows do begin
    for j=0, Cols do begin
      x=cos(Pi*((i-(Cols/2.0))*(1.0/Cols)))*cos(Pi*((j-(Rows/2.0))*(1.0/Rows)))
      hpMatrix[i*(Cols+1)+j]=(1.0-x)*(2.0-x)
    endfor
  endfor
  for i=0, Rows-1 do begin
    for j=0, Cols-1 do begin
      ARR[j,i]=temporary(ARR[j,i])*temporary(hpMatrix[i*(Cols+1)+j])
    endfor
  endfor
END



PRO ShiftArr, ARR
  ArrSize = SIZE(ARR)
  Cols = ArrSize[1]
  Rows = ArrSize[2]
  for i = 0, Rows - 1 do begin
    for j = 0, Cols - 1 do begin
      if ((i+j)mod 2) EQ 1 then begin
        Arr[j,i]=temporary(Arr[j,i]*(-1))
      endif
    endfor
  endfor
END


FUNCTION Normalize, ARR1
  ArrSize = SIZE(ARR1)
  Cols = ArrSize[1]
  Rows = ArrSize[2]
  ARR2 = temporary (FLTARR(COLS, ROWS,/Nozero))
  ARR2 = temporary(ARR1/255.0)
  RETURN, ARR2
END



PRO SRS_ENVIdoRegist_Start,event
;event handler for "doRegist" button:
forward_function ENVI_GET_DATA
widget_control,event.top,get_uvalue=state
if (state.Bimfid EQ -1) or (state.imagefid EQ -1) or (state.outfile EQ '' and state.saveflag EQ 0B)  then begin
  dummy=dialog_message(['Baseimage and Shiftimage must be imported, and You have to enter an output','filename before doing registe'],/error)
  return
endif else if state.Bimfid NE -1 and state.imagefid NE -1 then begin

rep_str=['Computing the shifts,rotation and scaling']
ENVI_REPORT_INIT,rep_str,base=rep_base,/interupt, title='Do registe'
ENVI_REPORT_INC,rep_base,state.Bimnl

Bimps = FLTARR(state.Bimns, state.Bimnl,/Nozero)
imageps = FLTARR(state.imagens, state.imagenl,/Nozero)

Bimps=ENVI_GET_DATA(fid=state.Bimfid,dims=state.Bimdims,pos=state.Bimpos)
imageps=ENVI_GET_DATA(fid=state.imagefid,dims=state.imagedims,pos=state.imagepos)

Im1 = FLTARR(state.Bimns, state.Bimnl,/Nozero)
Im2 = FLTARR(state.imagens, state.imagenl,/Nozero)

Im1 = temporary(Normalize(Bimps))
Im2 = temporary(Normalize(imageps))

ShiftArr, Im1
ShiftArr, Im2

fft_im1 = FLTARR(state.Bimns, state.Bimnl,/Nozero)
fft_im2 = FLTARR(state.imagens, state.imagenl,/Nozero)
fft_im1=temporary(FFT(im1,-1))
fft_im2=temporary(FFT(im2,-1))


absF_im1 = FLTARR(state.Bimns, state.Bimnl,/Nozero)
absF_im2 = FLTARR(state.imagens, state.imagenl,/Nozero)
absF_im1=temporary(abs(fft_im1))
absF_im2=temporary(abs(fft_im2))

HighPass, absF_im1
HighPass, absF_im2

pIm1 = FLTARR(state.Bimns, state.Bimnl,/Nozero)
pIm2 = FLTARR(state.imagens, state.imagenl,/Nozero)
pIm1 = temporary(LogPolar(absF_im1))
pIm2 = temporary(LogPolar(absF_im2))

fft_polar1 = FLTARR(state.Bimns, state.Bimnl,/Nozero)
fft_polar2 = FLTARR(state.imagens, state.imagenl,/Nozero)
fft_polar1=temporary(fft(pIm1,-1))
fft_polar2=temporary(fft(pIm2,-1))

R = temporary(Ritio(fft_polar1, fft_polar2))
inv_R = temporary(FFT(R, 1))

abs_IR = temporary(abs(inv_R))

maxn = MAX(abs_IR, position)
print,"maxn="
print, maxn
ArrSize = SIZE(abs_IR)
cols = ArrSize[1]
rows = ArrSize[2]
num = ArrSize[4]

b = 10 ^ (alog10(Cols) / Cols)
state.scale = b^(position MOD rows)
state.angle = 180.0 * (position / rows) / cols
ang=state.angle
if (ang eq 0.0) then begin
  widget_control, state.angtext, set_value=strtrim(string(state.angle),1)
  widget_control, state.scaltext, set_value=strtrim(string(state.scale),1)
   endif else if (ang gt 90.0) then begin
    ang=state.angle+180.0
     widget_control, state.angtext, set_value=strtrim(string(180.0-state.angle),1)
     widget_control, state.scaltext, set_value=strtrim(string(state.scale),1)
      endif else if (ang lt 90.0)  then begin
        ang=state.angle
        widget_control, state.angtext, set_value=strtrim(string(-state.angle),1)
        widget_control, state.scaltext, set_value=strtrim(string(state.scale),1)
      endif



Im3 = FLTARR(state.imagens, state.imagenl,/Nozero)
Im3=temporary(ROT(imageps, -ang, 1.0/state.scale, /CUBIC))

Ffft_im1=temporary(fft(Bimps,-1))
Ffft_im2=temporary(fft(Im3,-1))


R1=(Ffft_im1*temporary(conj(Ffft_im2)))
R2=(temporary(abs(Ffft_im1))*temporary(abs(Ffft_im2)))
R=temporary(R1)/temporary(R2)

IR=temporary(fft(R, 1))


maxn=MAX(IR,I)
print,"maxn="
print, maxn
state.IX=I MOD state.imagens
state.IY=I / state.imagens

X=state.IX
Y=state.IY

if (X gt 0.5*state.imagens) and ( Y gt 0.5*state.imagenl) then begin
   X=X-state.imagens
   Y=Y-state.imagenl
   widget_control, state.xtext, set_value=strtrim(string(X),1)
   widget_control, state.ytext, set_value=strtrim(string(Y),1)
   ;Im3=temporary (shift(Im3, X, Y))
endif else if (X gt 0.5*state.imagens) then begin
   X=X-state.imagens
   widget_control, state.xtext, set_value=strtrim(string(X),1)
   widget_control, state.ytext, set_value=strtrim(string(state.IY),1)
   ;Im3=temporary (shift(Im3, X, Y))
endif else if ( Y gt 0.5*state.imagenl) then begin
   Y=Y-state.imagenl
   widget_control, state.xtext, set_value=strtrim(string(state.IX),1)
   widget_control, state.ytext, set_value=strtrim(string(Y),1)
   ;Im3=temporary (shift(Im3, X, Y))
   endif else begin
     widget_control, state.xtext, set_value=strtrim(string(state.IX),1)
     widget_control, state.ytext, set_value=strtrim(string(state.IY),1)
    ; Im3=temporary (shift(Im3, X, Y))
 endelse
endif

Im3=temporary (shift(Im3, X, Y))

ENVI_REPORT_INIT,base=rep_base,/finish


rep_str=['Saving the warped image to the new file']
ENVI_REPORT_INIT,rep_str,base=rep_base,/interupt,title='warp image'
ENVI_REPORT_INC,rep_base,1

warpim=temporary(Im3)

;get data dimensions:
data_size=size(warpim,/dimensions)
ncols=data_size[0]
nrows=data_size[1]

;output the result:
bnames=['Warped image']
if state.saveflag EQ 0B then begin ;output to file:
  openw,lun,state.outfile,/GET_LUN
  writeu,lun,warpim
  close,lun
  free_lun,lun
  if state.mapflag EQ 1 then begin ;include georeferencing information:
    if state.psflag EQ 1 then envi_setup_head,bnames=bnames,data_type=1,$
      descrip='Warp result',file_type=0,interleave=0,nb=1,nl=nrows,$
      ns=ncols,/open,r_fid=warpim_fid,sensor_type=no,fname=state.outfile,/write,$
      map_info=state.imagemap,pixel_size=state.imageps ELSE envi_setup_head,bnames=bnames,data_type=1,$
      descrip='Warp result',file_type=0,interleave=0,nb=1,nl=nrows,$
      ns=ncols,/open,r_fid=warpim_fid,sensor_type=no,fname=state.outfile,/write,$
      map_info=state.Bimmap
    endif
  if state.mapflag EQ 0 then begin ;exclude georeferencing information:
    if state.psflag EQ 1 then envi_setup_head,bnames=bnames,data_type=1,$
      descrip='warp result',file_type=0,interleave=0,nb=1,nl=nrows,$
      ns=ncols,/open,r_fid=rgb_fid,sensor_type=no,fname=state.outfile,/write,$
      pixel_size=state.imageps ELSE envi_setup_head,bnames=bnames,data_type=1,$
      descrip='Warp result',file_type=0,interleave=0,nb=1,nl=nrows,$
      ns=ncols,/open,r_fid=warpim_fid,sensor_type=no,fname=state.outfile,/write,$
      def_stretch={dstretch_struct,type:1L,vals:[2.0,0.0,2.0,98.0]}
    endif
endif else begin ;output to memory:
  if state.mapflag EQ 1 then envi_enter_data,warpim,$
    descrip='Warp result',sensor_type=no,bnames=bnames,$
    file_type=0,map_info=state.Bimmap $
  ELSE envi_enter_data,warpim,$
    descrip='Warp result',sensor_type=no,bnames=bnames,$
    file_type=0
endelse
;kill status report dialog:
ENVI_REPORT_INIT,base=rep_base,/finish
END


PRO SRS_ENVIdoRegist_Cancel,event
;event handler for "Cancel" button:
widget_control,event.top,get_uvalue=state
device,decomposed=state.olddc
widget_control,event.top,/destroy
return
END


PRO SRS_ENVIRegisterLAbout,event
;event handler for "About this User Function" button:
dummy=dialog_message(['SRS_ENVI (Shift, Rotation, Scaling) written by IDL/ENVI Team at UTEP.',$
                      'It is used for two images that shift, rotation and scaling each other,',$
                      'but the overlapping parts must be larger than 30%. the biggest scaling ',$
                      'times can be 1.85, the rotation angle can be any degree, but for simplicity,',$
                      'we usually set the angle in (-90, 90). If your base image is geocoded,',$
                      'the map information will automatically apply to the registered (warped) image',$
                      '                                  ',$
                      'Several Comments about how to use:',$
                      '       1.  The size of image can be any dimensions,but should',$
                      '           be exactly the same of two images and square,',$
                      '       2.  For avoiding the confuse, it is better to',$
                      '           set two images Upper Left to (1,1)',$                      '
                      '       3.  After hit the doRegiste button, the rotation, ',$
                      '           scaling times and shifts will be shown on',$
                      '           a new image with no rotation, no scaling, no shifts',$
                      '           can be gotten',$
                      '       4.  The exact shift value (column, row),angle and scale will be used',$
                      '           while the user function do the actual registration.',$
                      '                                          ',$
                      'For more information about the algorithm or how to use',$
                      'Please contact IDL/ENVI Team at UTEP, or hongjiexie@yahoo.com',$
                      '      ',$
                      'Instructor:Drs. Kreinovich and Keller'],/info)

END


PRO NULL_EVENTS,event
;event handler for null events
END



PRO SRS_ENVI,event
;
;++++++++++++++++++++++
;MAIN ROUTINE PROCEDURE
;++++++++++++++++++++++
;
;specify ENVI functions so they are not treated as arrays being subscripted:
forward_function widget_outfm
;error checking:
!error_state.code=0
catch,error
if error NE 0 then begin
  help,/last_message,output=traceback
  errarray=['Error Caught',traceback]
  dummy=dialog_message(errarray,/error,/cancel)
  if STRUPCASE(dummy) EQ 'CANCEL' then return
endif

;main widget creation:
tlb=widget_base(title='AutoRegiste Made By IDL/ENVI Team at UTEP',/column,xoffset=125,yoffset=25,/tlb_kill_request_events)

group1=widget_base(tlb,/column,/frame)
  row1=widget_base(group1,/row)
    label1=widget_label(row1,VALUE='Select Base Image:')
  row2=widget_base(group1,/row)
    text1=widget_text(row2,xsize=40,ysize=1,EVENT_PRO='NULL_EVENTS')
  row3=widget_base(group1,/row)
    importbut1=widget_button(row3,VALUE='Import Band...',EVENT_PRO='SRS_ENVIOpenBaseim')
    deletebut1=widget_button(row3,VALUE='Remove',EVENT_PRO='SRS_ENVIRemoveBaseim')

group2=widget_base(tlb,/column,/frame)
  row4=widget_base(group2,/row)
    label2=widget_label(row4,VALUE='Select Image with Shift, Rotation and Scaling:')
  row5=widget_base(group2,/row)
    text2=widget_text(row5,xsize=40,ysize=1,EVENT_PRO='NULL_EVENTS')
  row6=widget_base(group2,/row)
    importbut2=widget_button(row6,VALUE='Import Band...',EVENT_PRO='OpenSRS_ENVIim')
    deletebut2=widget_button(row6,VALUE='Remove',EVENT_PRO='RemoveSRS_ENVIim')

group3=widget_base(tlb, /column,/frame)
  row7=widget_base(group3,/row)
     label3=widget_label(row7, VALUE='The computation results are: ')
  row8=widget_base(group3,/row)
     label4=widget_label(row8, VALUE='Rotation angle (degrees) = ')
     angtext=widget_text(row8, xsize=10, ysize=1, EVENT_PRO='NULL_EVENTS')
     label5=widget_label(row8, VALUE='Scaling (times) = ')
     scaltext=widget_text(row8, xsize=10, ysize=1, EVENT_PRO='NULL_EVENTS')
  row9=widget_base(group3,/row)
     label6=widget_label(row9, VALUE='Shift (pixels): Column = ')
     xtext=widget_text(row9, xsize=10, ysize=1, EVENT_PRO='NULL_EVENTS')
     label7=widget_label(row9, VALUE='Row = ')
     ytext=widget_text(row9, xsize=10, ysize=1, EVENT_PRO='NULL_EVENTS')
  row10=widget_base(group3,/row)
    fileio=widget_outfm(row10)
  row11=widget_base(group3,/row)
    start2but=widget_button(row11,VALUE='doRegist',EVENT_PRO='SRS_ENVIdoRegist_Start')
    cancel2but=widget_button(row11,VALUE='Cancel',EVENT_PRO='SRS_ENVIdoRegist_Cancel')
    aboutbut=widget_button(row11,VALUE='About this User Function',EVENT_PRO='SRS_ENVIRegisterLAbout')

widget_control,tlb,/realize
;acquire current color decomposition state:
device,get_decomposed=olddc
;adjust to 8-bit mode if necessary:
if !D.N_COLORS GT 256 then device,decomposed=0


;setup state structure:

state={text1:text1,IX:float(1), IY:float(1),angle:float(1), scale:float(1),xtext:xtext,ytext:ytext,angtext:angtext,scaltext:scaltext,Bimname:'',Bimpos:lonarr(1),Bimfid:-1L,Bimns:0L,Bimnl:0L,$
       Bimps:[0.0,0.0],imageps:[0.0,0.0],text2:text2,imagefid:-1L,imagepos:lonarr(1),imagens:0L,$
       imagenl:0L,elev:30,azim:0,tlb:tlb,tlb2:0L,elevtext:0L,azimtext:0L,xsizetext:0L,dembs:0,$
       imagedims:lonarr(5),ysizetext:0L,xddtext:0L,yddtext:0L,olddc:olddc,previewflag:0B,$
       Bimdims:lonarr(5),draw:0L,outfile:'',imagename:'',fileio:fileio,$
       metext:0L,saveflag:0B,preview_id:0L,Bim_ptr:ptr_new(/allocate_heap),tlb5:0L,$
       Bimmap:{envi_map_struct},imagemap:{envi_map_struct}, stretchflag:0B,mapflag:0B,psflag:0B}
widget_control,tlb,set_uvalue=state
;register widget with XMANAGER:
XMANAGER,'SRS_ENVI',tlb ; , /NO_BLOCK


END
