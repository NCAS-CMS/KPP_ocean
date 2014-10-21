SUBROUTINE mckpp_output_control(kpp_3d_fields,kpp_const_fields,kpp_timer)
  IMPLICIT NONE

#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  TYPE(kpp_timer_type) :: kpp_timer

  ! Need to drive MEAN_OUTPUT and RANGE_OUTPUT from ncdf_out.f
  ! Need to handle rotation of all output files
  ! Need to handle writing of all output

END SUBROUTINE mckpp_output_control

SUBROUTINE mckpp_output_open(file,ncid)      
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
  !#include <kpp_3d_type.com>
  !#include <output.com>

  INTEGER status
  INTEGER,intent(out) :: ncid
  CHARACTER*50,intent(in) :: file      
  
  status=NF_OPEN(file,NF_WRITE,ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  RETURN
END SUBROUTINE mckpp_output_open

SUBROUTINE mckpp_output_close(ncid)      
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
  !#include <kpp_3d_type.com>
  !#include <output.com>

  INTEGER status
  INTEGER, intent(in) :: ncid
  
  status=NF_CLOSE(ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  RETURN
END SUBROUTINE mckpp_output_close
