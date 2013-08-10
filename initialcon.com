c      REAL hu(NPTS),SSU(NPTS,NVEL),hx(NPTS),SSX(NPTS,NSCLR)
c      INTEGER mu(NPTS),mx(NPTS)
c      CHARACTER*40 initprofile
c      common/ initialcon/ hu,SSU,hx,mu,mx,SSX,
c     +                    initprofile

      LOGICAL L_INITDATA,L_INTERPINIT,L_RESTART
      CHARACTER*40 initdata_file,restart_infile
      common/ initdata / L_INITDATA,L_INTERPINIT,L_RESTART,
     +                  initdata_file,restart_infile
      
