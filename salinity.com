      LOGICAL L_RELAX_SAL, L_UPD_SAL, L_PERIODIC_SAL
      INTEGER ndtupdsal, sal_period
      REAL sal_clim(NPTS,NZP1),scorr(NPTS,NZP1),relax_sal_in(NY),
     +     relax_sal(NPTS)
      CHARACTER*40 sal_file
      COMMON /sal_in/ L_RELAX_SAL, L_UPD_SAL, ndtupdsal, sal_file,
     +     L_PERIODIC_SAL, sal_period, scorr, relax_sal_in, relax_sal
