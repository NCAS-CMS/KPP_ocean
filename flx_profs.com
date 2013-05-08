      REAL wU(NPTS,0:NZtmax,NVP1),wX(NPTS,0:NZtmax,NSP1),
     +     wXNT(NPTS,0:NZtmax,NSCLR),
     +     wtI(NPTS,0:NZ,0:NJDT),wsI(NPTS,0:NZ,0:NJDT) 
      common/ flx profs / wU,wX,
     +     wXNT,
     +     wtI,wsI

