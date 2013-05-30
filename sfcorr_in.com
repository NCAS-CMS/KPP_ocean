      LOGICAL L_SFCORR_WITHZ, L_SFCORR, L_UPD_SFCORR, L_PERIODIC_SFCORR 
      INTEGER ndtupdsfcorr, sfcorr_period
      REAL sfcorr_withz(NPTS,NZP1), Sinc_fcorr(NPTS,NZP1),
     &     sfcorr_twod(NPTS)
      CHARACTER*40 sfcorrin_file
      common /sfcorr_in/ L_SFCORR_WITHZ,L_SFCORR,
     &     ndtupdsfcorr, sfcorrin_file,L_UPD_SFCORR, L_PERIODIC_SFCORR,
     &     sfcorr_period, sfcorr_withz, Sinc_fcorr
