SUBROUTINE MCKPP_ABORT

  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#ifdef COUPLE
#ifdef OASIS2
  ! Support for stopping in OASIS2
  ! NPK April 2009 - R2
  CALL halte('KPP : MIXED_ABORT called')
#else
  ! Support for stopping in OASIS3
  ! NPK 2/11/09 - R3
#ifdef OASIS3
  CALL mpi1_oasis3_terminate()
#else
#ifdef CFS
  ! Support for stopping with the CFS coupler
  ! Unsure how to stop the model for the GFS - Just stop?
  ! NPK June 2009 - R2     
  STOP
#endif /*CFS*/
#endif /*OASIS3*/
#endif /*OASIS2*/
#else
  STOP
#endif /*COUPLE*/
  
END SUBROUTINE MCKPP_ABORT
