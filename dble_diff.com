      REAL dift(npts,0:NZtmax),Rp(NZtmax)
      common/ dble diff / dift,Rp

