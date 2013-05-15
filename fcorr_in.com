      LOGICAL L_FCORR_WITHZ, L_FCORR, L_UPD_FCORR, L_PERIODIC_FCORR, 
     &     L_NO_FREEZE, L_NO_ISOTHERM
      INTEGER ndtupdfcorr, fcorr_period, isotherm_bottom
      REAL fcorr_withz(NPTS,NZP1), Tinc_fcorr(NPTS,NZP1),
     &     fcorr_twod(NPTS), isotherm_threshold
      CHARACTER*40 fcorrin_file
      common /fcorr_in/ L_FCORR_WITHZ,L_FCORR,
     &     ndtupdfcorr, fcorrin_file,L_UPD_FCORR, L_PERIODIC_FCORR,
     &     fcorr_period, Tinc_fcorr, fcorr_withz, L_NO_FREEZE,
     &     L_NO_ISOTHERM, isotherm_bottom, isotherm_threshold
