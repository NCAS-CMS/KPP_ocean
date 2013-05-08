      INTEGER kmix(npts)
      REAL hmix(npts),difm(npts,0:NZtmax),
     $          difs(npts,0:NZtmax),ghat(npts,NZtmax)
      common/ kprof out / hmix,difm,difs,ghat,kmix

