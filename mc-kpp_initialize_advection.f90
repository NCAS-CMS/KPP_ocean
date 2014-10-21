SUBROUTINE mckpp_initialize_advection(kpp_3d_fields)
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
#include <ocn_advec.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  INTEGER nmode(npts)
  
  INTEGER ipt,ivar,imode,status
  
  IF (L_ADVECT) THEN
     status=NF_OPEN(advect_file,0,ncid_advec)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)     
     call MCKPP_READ_IPAR(kpp_3d_fields,ncid_advec,'nmode_tadv',1,1,kpp_3d_fields%nmodeadv(:,1))     
     call MCKPP_READ_IPAR(kpp_3d_fields,ncid_advec,'mode_tadv',maxmodeadv,1,kpp_3d_fields%modeadv(:,:,1))
     call MCKPP_READ_PAR(kpp_3d_fields,ncid_advec,'tadv',maxmodeadv,1,kpp_3d_fields%advection(:,:,1))
     call MCKPP_READ_IPAR(kpp_3d_fields,ncid_advec,'nmode_sadv',1,1,kpp_3d_fields%nmodeadv(:,2))
     call MCKPP_READ_IPAR(kpp_3d_fields,ncid_advec,'mode_sadv',maxmodeadv,1,kpp_3d_fields%modeadv(:,:,2))
     call MCKPP_READ_PAR(kpp_3d_fields,ncid_advec,'sadv',maxmodeadv,1,kpp_3d_fields%advection(:,:,2))  
     status=NF_CLOSE(ncid_advec)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ELSE
     DO ipt=1,npts
        DO ivar=1,2
           kpp_3d_fields%nmodeadv(ipt,ivar)=0
        ENDDO
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE mckpp_initialize_advection
