      REAL CP(NPTS,0:NZP1tmax),rho(NPTS,0:NZP1tmax),rhoh2o(NPTS),
     &     rhob,talpha(0:NZP1tmax),sbeta(0:NZP1tmax),
     &     Sref(NPTS),SSref(NPTS),epsw  
      LOGICAL L_SSRef
      common/ ocn_paras / CP,rho,rhoh2o,rhob,   
     &     talpha,sbeta,Sref,SSref,epsw,L_SSref

