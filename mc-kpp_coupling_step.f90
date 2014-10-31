#ifdef COUPLE
SUBROUTINE mckpp_coupling_step(kpp_3d_fields,kpp_const_fields)
  
  ! Wrapper to choose correct method to return oceanic fields to atmosphere via coupler
  ! Called only #ifdef COUPLE
  
  IMPLICIT NONE
#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
#ifdef OASIS3
  CALL MCKPP_COUPLING_OASIS3_OUTPUT(kpp_3d_fields,kpp_const_fields)
#endif /*OASIS3*/
  
  RETURN
END SUBROUTINE mckpp_coupling_step
#endif /*COUPLE*/
