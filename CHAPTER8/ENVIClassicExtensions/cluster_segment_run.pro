; docformat = 'rst'
; cluster_segment_run.pro
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

PRO cluster_segment_run_define_buttons, buttonInfo
   ENVI_DEFINE_MENU_BUTTON, buttonInfo, $
      VALUE = 'Cluster Segments', $
      REF_VALUE = 'K-Means', $
      EVENT_PRO = 'cluster_segment_run', $
      UVALUE = 'SEGMENT_CLASS',$
      POSITION = 'after'
END

;+
; :Description:
;       ENVI extension for unsupervised classification
;       of a segmented image using Hu invariant moments
;       with an agglomerative hierarchical clustering algorithm 
; :Params:
;      event:  in, required 
;         if called from the ENVI menu                  
; :Uses:
;       ENVI
;       COYOTE
;       HU_MOMENTS
;       HCL 
;       CLASS_LOOKUP_TABLE
;       SEGMENT_CLASS
; :Author:
;      Mort Canty (2013)       
;-
Pro cluster_segment_run, event

COMPILE_OPT IDL2

   print, '---------------------------------------------'
   print, ' Segment Classification with Hu moments'
   print, systime(0)
   print, '---------------------------------------------'
   
   catch, theError
   if theError ne 0 then begin
      void = Dialog_Message(!Error_State.Msg, /error)
      return
   endif   

; select segment or classification file
   envi_select,fid=fid,pos=pos,dims=dims,/no_dims,/no_spec,/band_only,title='Select segment image'
   if fid eq -1 then begin
      print, 'cancelled'
      return
   endif
   
   envi_file_query, fid, file_type=file_type,class_names=class_names,$
             data_type=data_type,fname=fname,xstart=xstart,ystart=ystart
   
   num_cols = dims[2]-dims[1]+1
   num_rows = dims[4]-dims[3]+1
   
   climg = envi_get_data(fid=fid,dims=dims,pos=pos)
   
   if file_type eq 3 then begin
; this is a classification file   
      base = widget_auto_base(title='Select classes')  
      list = class_names[1:*] 
      wm = widget_multi(base, list=list, uvalue='list', /auto)  
      result = auto_wid_mng(base)  
      if (result.accept eq 0) then begin
         print, 'cancelled'
         return
      endif  
;    the K selected classes
      cl = where(result.list eq 1,K)+1
   end else begin
; this is a segemented file   
      K = max(climg)
      cl = lindgen(K)+1
   endelse      
   
; tie point
   map_info = envi_get_map_info(fid=fid)
   envi_convert_file_coordinates, fid, dims[1], dims[3], e, n, /to_map
   map_info.mc[2:3]= [e,n]   
   
; number of clusters
   base = widget_auto_base(title='Number of Classes')
   wg = widget_sslider(base, title='Classes', min=2, max=15, $
     value=5, dt=1, uvalue='slide', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
      print, 'cancelled'
      return
   endif
   K = fix(result.slide)
   print, 'Number of Classes', K   
   
; minimum segment size
   base = widget_auto_base(title='Minimum size')
   wg = widget_sslider(base, title='SIze', min=0, max=500, $
     value=100, dt=2, uvalue='slide', /auto)
   result = auto_wid_mng(base)
   if (result.accept eq 0) then begin
      print, 'cancelled'
      return
   endif
   min_size = result.slide   
   
; blobbing for connected segments   
   segmented = segment_class(climg,cl,min_size=min_size,/all_neighbors)   
   num_segs = max(segmented)
      
; calculate table of invariant moments
   moments = fltarr(7,num_segs)
   progressbar = Obj_New('cgprogressbar',/cancel,$
                    title='Calculating moments...')
   progressbar->start
   window, 11,xsize=200,ysize=200,title='Segment'
   A0=fltarr(5,5)
   for i=0L, num_segs-1 do begin
      if progressbar->CheckCancel() then begin
         progressbar->destroy
         wdelete,11
         return
      endif
      pct = (i)*100.0/(num_segs)
      progressbar->Update,fix(pct)
      indices = where(segmented eq i+1,count)
      if count gt 0 then begin
         X = indices mod num_cols
         Y = indices/num_cols
         A = segmented[min(X):max(X)+1,min(Y):max(Y)+1]
         indices1 = where(A eq i+1,complement=complement,ncomplement=ncomplement)
         A[indices1]=1.0 
         if ncomplement gt 0 then A[complement]=0.0
         if i mod 5 eq 0 then begin
             wset,11
             tvscl,A0
             tvscl,A
             A0=A*0
         endif
         moments[*,i] = hu_moments(A,/log)         
      endif
   endfor
   progressbar->destroy
   wdelete,11
   
   moments = standardize(temporary(moments)) 
        
; cluster
   HCL,moments,K,Ls
   if n_elements(Ls) eq 1 then return   
   
; plot in 3D feature space   
   symbols = objarr(num_segs)
   for i=0L,num_segs-1 do begin
      color=class_lookup_table(Ls[i]+1) 
      oOrb = OBJ_NEW('orb', color=color)
      oOrb->Scale, 0.1, 0.1, 0.1 
      symbols[i] = OBJ_NEW('IDLgrSymbol', oOrb)
   endfor    
   moments = (moments = moments > (-5)) < 5 
   xplot3d,moments[0,*],moments[1,*],moments[2,*],linestyle=6,symbol=symbols,$
      xrange=[-5,5],yrange=[-5,5],zrange=[-5,5], $
      title='Hu Moments', xtitle='log(h1)',ytitle='log(h2)',ztitle='log(h3)'   
   
; re-label segments with cluster label
   _ = histogram(segmented,reverse_indices=r)
   segmented=segmented*0
   for i=1L,num_segs do begin
       p = r[r[i]:r[i+1]-1]   ;pixels in segment 
       segmented[p]=Ls[i-1]+1 ;new label
   endfor
   
   c_names=strarr(K+1)
   c_names[0]='unclassified'
   for i=1,K do c_names[i]='cluster: '+strtrim(i,2)
   
   envi_enter_data, segmented, $
      file_type=3, $
      map_info=map_info, $
      xstart=xstart+dims[1], $
      ystart=ystart+dims[3], $
      bnames=['EM-SEG('+fname+')'],$
      num_classes=K+1, $
      class_names=c_names, $
      lookup=class_lookup_table(indgen(K+1))

End