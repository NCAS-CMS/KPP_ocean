      LOGICAL LKPP , LRI, LDD, LICE, LBIO, LNBFLX, LTGRID, LRHS,
     &     lsaveaverages,L_STRETCHGRID,lrepeatdat,
     &     lradtq,lfluxSSTdat,
     &     lLWupSSTdat,lrhdat,lclddat,L_EKMAN_PUMP,
     &     L_BARRIER_REMOVE,L_BARRIER_SALISO,L_BARRIER_SALVAVG,
     &     L_NO_EGTP
  
      INTEGER barrier_subdepth,barrier_dT,barrier_ifirst,barrier_ilast,
     &     barrier_jfirst,barrier_jlast
      
      common/ proc_swit / LKPP , LRI, LDD, LICE, LBIO, LNBFLX, LTGRID,
     &     LRHS, 
     &     lsaveaverages,L_STRETCHGRID,lrepeatdat,
     &     lradtq,lfluxSSTdat,
     &     lLWupSSTdat,lrhdat,lclddat,L_EKMAN_PUMP

