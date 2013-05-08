      REAL hu(NPTS),SSU(NPTS,NVEL),hx(NPTS),SSX(NPTS,NSCLR)
      INTEGER mu(NPTS),mx(NPTS)
      CHARACTER*40 initprofile,restart_infile
      common/ initialcon/ hu,SSU,hx,mu,mx,SSX,
     +                    initprofile

      LOGICAL L_INITDATA,L_INTERPINIT,L_RESTART
      CHARACTER*40 initdata_file
      common/ initdata / L_INITDATA,L_INTERPINIT,L_RESTART,
     +                  initdata_file,restart_infile
      
