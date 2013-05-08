      LOGICAL L_RELAX_SAL, L_UPD_SAL, L_PERIODIC_SAL
      INTEGER ndtupdsal, sal_period
      REAL relax_sal_in(NY),relax_sal(NPTS)
      CHARACTER*40 sal_file
      COMMON /sal_in/ L_RELAX_SAL, L_UPD_SAL, ndtupdsal, sal_file,
     +     L_PERIODIC_SAL, sal_period, relax_sal_in, relax_sal

      LOGICAL L_RELAX_OCNT,L_UPD_OCNT,L_PERIODIC_OCNT
      INTEGER ndtupdocnt, ocnt_period
      REAL relax_ocnT_in(NY),relax_ocnT(NPTS)
      CHARACTER*40 ocnT_file
      COMMON /ocnT_in/ L_RELAX_OCNT,L_UPD_OCNT,ndtupdocnt,ocnt_period,
     +     L_PERIODIC_OCNT,ocnt_file,relax_ocnT_in,relax_ocnT
