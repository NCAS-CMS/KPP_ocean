      REAL DMAX
      LOGICAL L_VGRID_FILE
      CHARACTER(LEN=50) vgrid_file
      common/ vert pgrid/ DMAX,
     +     L_VGRID_FILE,vgrid_file

