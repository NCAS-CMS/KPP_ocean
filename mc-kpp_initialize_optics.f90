SUBROUTINE mckpp_initialize_optics(kpp_3d_fields)
  
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)

#include <netcdf.inc>  
#include <mc-kpp_3d_type.com>
#include <proc_pars.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  integer jerlov(npts)
  integer ipt, status  
  LOGICAL L_PARAS
  
  L_PARAS=L_JERLOV

  IF (L_PARAS) THEN
     status=NF_OPEN(paras_file,0,ncid_paras)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ENDIF
  
  IF (L_JERLOV) THEN
     call MCKPP_READ_IPAR(kpp_3d_fields,ncid_paras,'jerlov',1,1,kpp_3d_fields%jerlov)
  ELSE
     DO ipt=1,npts
        kpp_3d_fields%jerlov(ipt)=3
     ENDDO
  ENDIF
  
  IF (L_PARAS) THEN
     status=NF_CLOSE(ncid_paras)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ENDIF
  
  RETURN
END SUBROUTINE mckpp_initialize_optics
