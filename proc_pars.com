c      INTEGER jerlov(npts)
c      common/ proc_pars / jerlov
      
      CHARACTER*40 paras_file
      INTEGER ncid_paras,max_ekman_depth,max_ekadv_depth
      LOGICAL L_JERLOV
      common/ l_proc_pars / paras_file,ncid_paras,
     $     l_jerlov,max_ekman_depth,max_ekadv_depth
      
