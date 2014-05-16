; $Id: //depot/Release/ENVI50_IDL82/idl/idldir/lib/clust_wts.pro#1 $
;
; Copyright (c) 1996-2012, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
; NAME:
;       CLUST_WTS1
;
; PURPOSE:
;       This function computes the weights (the cluster centers) of an
;       M-column, N-row array, where M is the number of variables and
;       N is the number of observations or samples. The result is an
;       M-column, N_CLUSTERS-row array of cluster centers.
;       
;       mModified to include CGPROGRESSBAR my Mort Canty (2013)
;
; CATEGORY:
;       Statistics
;
; CALLING SEQUENCE:
;       Result = Clust_wts(Array)
;
; INPUTS:
;       Array:    An M-column, N-row array of type float or double.
;
; KEYWORD PARAMETERS:
;             DOUBLE:  If set to a non-zero value, computations are done in
;                      double precision arithmetic.
;
;         N_CLUSTERS:  Use this keyword to specify the number of cluster
;                      centers. The default is to compute N cluster centers.
;
;       N_ITERATIONS:  Use this keyword to specify the number of iterations
;                      in computing the cluster centers. The default is to
;                      use 20 iterations.
;
;       VARIABLE_WTS:  An M-element vector of variable weights. The elements
;                      of this vector are used to give greater or lesser
;                      importance to each variable (each column) in determining
;                      the cluster centers. The default is to give all variables
;                      equal weighting using a value of 1.0.
;
; EXAMPLE:  See the documentation for CLUSTER.
;
; REFERENCE:
;       CLUSTER ANALYSIS (Third Edition)
;       Brian S. Everitt
;       ISBN 0-340-58479-3
;
; MODIFICATION HISTORY:
;       Written by:  GGS, RSI, June 1996
;                    Adapted from an algorithm written by Robb Habbersett
;                    of Los Alamos National Laboratory.
;-

FUNCTION Clust_Wts1, Array, Double = Double, N_Clusters = N_Clusters, $
                           N_Iterations = N_Iterations, Variable_Wts = Variable_Wts

  ON_ERROR, 2

  Dimension = SIZE(Array)
  if Dimension[0] ne 2 then MESSAGE, "Input array must be a two-dimensional."
  N_Variables = Dimension[1]
  N_Samples = Dimension[2]

  ;Set default keyword values.
  if N_ELEMENTS(Double) eq 0 then Double = (Dimension[Dimension[0]+1] eq 5)
  Zero = Double ? 0.0d : 0.0

  if N_ELEMENTS(N_Clusters)   eq 0 then N_Clusters   = Dimension[2]
  if N_ELEMENTS(N_Iterations) eq 0 then N_Iterations = 20

  ;Weighting of variables.
  if KEYWORD_SET(Variable_Wts) eq 0 then $
    Variable_Wts = REPLICATE(1.0+Zero, N_Variables) $
  else if N_ELEMENTS(Variable_Wts) ne N_Variables then $
    MESSAGE, "N_Variables parameter must be an M-element vector."

  ;Initial 'learning rate'.
  Learning = [0.5, 0.1] + Zero
  InitialLR = Learning[0]

  DeltaLR = (Learning[0] - Learning[1]) / N_Iterations

  ;Work arrays.
  WorkCol = REPLICATE(1.0 + Zero, 1, N_Clusters)
  WorkRow = REPLICATE(1.0 + Zero, 1, N_Variables)

  Count = LONARR(N_Clusters)

  ;Normalized uniformly random cluster weights.
  Weights = RANDOMU(SEED, N_Variables, N_Clusters) + Zero
  for k = 0L, N_Clusters-1 do Weights[*,k] = $
    (Weights[*,k] / TOTAL(Weights[*,k])) * Variable_Wts
    
   progressbar = Obj_New('cgprogressbar',$
             title='Clustering ... ',xsize=300,ysize=20,/cancel)            
   progressbar->start    

  for k = 1L, N_Iterations do begin
    if progressbar->CheckCancel() then begin
       print,'Clustering aborting ...'
       k = N_Iterations
    endif
    pct=(k-1)*100/N_Iterations      
    progressbar->update,fix(pct)             
    for Sample = 0L, N_Samples-1 do begin
      Vector = Array[*,Sample] # WorkCol - Weights
      Metric = WorkRow # ABS(Vector)
      Minimal = WHERE(Metric eq MIN(Metric))
      if n_elements(minimal) eq 1 then begin
          ; Need to convert to a scalar so we can use our "0" index
          ; trick below.
          Minimal = Minimal[0]
          Count[Minimal] = Count[Minimal] + 1 ;Optional output parameter.
          ; We use a 0 to subscript into the Weights column index.
          ; IDL will automatically copy the entire row so this is the
          ; same as using a "*" but avoids the overhead of computing
          ; the * indices and speeds up the loop by ~20%.
          Weights[0,Minimal] = InitialLR * Vector[*,Minimal] + $
            Weights[*,Minimal]
      endif else begin
          Count[Minimal] = Count[Minimal] + 1 ;Optional output parameter.
          Weights[*,Minimal] = InitialLR * Vector[*,Minimal] + $
            Weights[*,Minimal]
      endelse

    endfor

    InitialLR = InitialLR - DeltaLR
  endfor
  progressbar->destroy
  RETURN, Weights
          ;This procedure initializes a potentially different set of weights
          ;each time it is called; Even if the input data is identical between
          ;callings.

END
