      REAL focn(npts),Tref(npts),Ssurf(npts),uref(npts),vref(npts),
     $    Qs(npts),ug0(npts),vg0(npts),albocn
      common/ ocn state / focn,Tref,Ssurf,uref,vref,Qs,ug0,vg0,albocn

