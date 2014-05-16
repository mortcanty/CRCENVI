; docformat = 'rst'
; cluster_fkm.pro
;+
; :Description:
;      Distance clusterer from IDL library.
;      Modified for FKM (Mort Canty (2006))      
;-
FUNCTION Cluster_fkm, Array, Weights, Double = Double, N_clusters = N_clusters
; modified distance clusterer from IDL library
  ON_ERROR, 2
  Dimension = SIZE(Array)
  if Dimension[0] ne 2 then MESSAGE, "Input array must be a two-dimensional."
  if N_ELEMENTS(Double) eq 0 then Double = (Dimension[Dimension[0]+1] eq 5)
  if Double eq 0 then Zero = 0.0 else Zero = 0.0d ;Type casting constant.
  if KEYWORD_SET(N_Clusters) eq 0 then N_Clusters = ((SIZE(Weights))[2])
  ;Work arrays.
  WorkRow = REPLICATE(1.0, 1, Dimension[1]) + Zero
  WorkCol = REPLICATE(1.0, 1, N_Clusters) + Zero
  ClusterNumber = LONARR(Dimension[2])
  progressbar = Obj_New('progressbar', Color='blue', Text='0',$
              title='classifying...',xsize=250,ysize=20)
  progressbar->start
  for Sample = 0L, Dimension[2]-1 do begin
    if (sample mod 1000) eq 0 then begin
       if progressbar->CheckCancel() then begin
          print,'aborted'
          progressbar->Destroy
          return, 0
       endif
       pct=sample*100/(Dimension[2]-1)
       progressbar->Update,pct,text=strtrim(pct,2)+'%'
    endif
    Vector = Array[*,Sample] # WorkCol - Weights
    Metric = WorkRow # ABS(Vector)
    ClusterNumber[Sample] = (WHERE(Metric eq MIN(Metric)))[0]
  endfor
  progressbar->destroy
  RETURN, TRANSPOSE(ClusterNumber)
END