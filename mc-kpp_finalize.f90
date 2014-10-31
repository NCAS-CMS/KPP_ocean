SUBROUTINE mckpp_finalize(kpp_3d_fields,kpp_const_fields,kpp_timer)
  IMPLICIT NONE

#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  TYPE(kpp_timer_type) :: kpp_timer

  ! Checkpoint if requested
  CALL MCKPP_RESTART_CONTROL(kpp_3d_fields,kpp_const_fields,kpp_timer)

  ! Close output files
  IF (kpp_const_fields%L_OUTPUT_INST) &
       CALL MCKPP_OUTPUT_CLOSE(kpp_const_fields%ncid_out)
  IF (kpp_const_fields%L_OUTPUT_MEAN) &
       CALL MCKPP_OUTPUT_CLOSE(kpp_const_fields%mean_ncid_out)
  IF (kpp_const_fields%L_OUTPUT_RANGE) THEN
     CALL MCKPP_OUTPUT_CLOSE(kpp_const_fields%min_ncid_out)
     CALL MCKPP_OUTPUT_CLOSE(kpp_const_fields%max_ncid_out)
  ENDIF
  
  ! Terminate coupled model if coupled
#ifdef COUPLE
#ifdef OASIS3
  CALL MCKPP_COUPLING_OASIS3_FINALIZE()
#endif /*OASIS3*/
#endif /*COUPLE*/
  
  !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
  !CALL KPP_TIMER_PRINT(kpp_timer)

  RETURN
END SUBROUTINE mckpp_finalize
