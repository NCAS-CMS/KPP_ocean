      REAL DMAX,zm(NZP1),hm(NZP1),dm(0:NZ)
      LOGICAL L_VGRID_FILE
      CHARACTER(LEN=50) vgrid_file
      common/ vert pgrid/ DMAX,zm,hm,dm,
     +     L_VGRID_FILE,vgrid_file

