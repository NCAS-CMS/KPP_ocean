      REAL dlat(npts),dlon(npts),rlat(npts),rlon(npts),f(npts)
      LOGICAL L_REGGRID
      common/ location  / dlat,dlon,rlat,rlon,f,L_REGGRID

