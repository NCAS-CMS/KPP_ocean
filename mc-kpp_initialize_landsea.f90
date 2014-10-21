SUBROUTINE mckpp_initialize_landsea(kpp_3d_fields)
  
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
 
#include <netcdf.inc>
! Automatically includes parameter.inc
#include <mc-kpp_3d_type.com>
#include <landsea.com>

  REAL :: landsea(npts)
  TYPE (kpp_3d_type) :: kpp_3d_fields
  
  INTEGER :: ipt,status
  
  IF (L_LANDSEA) THEN
     status=NF_OPEN(landsea_file,0,ncid_landsea)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)     
     WRITE(6,*) 'Initialized landsea file',ncid_landsea
     
     call MCKPP_READ_PAR(kpp_3d_fields,ncid_landsea,'lsm',1,1,landsea)
     DO ipt=1,npts
        IF (landsea(ipt) .EQ. 1.0) THEN
           kpp_3d_fields%L_OCEAN(ipt)=.FALSE.
        ELSE
           kpp_3d_fields%L_OCEAN(ipt)=.TRUE.
        ENDIF
     ENDDO
     
     call MCKPP_READ_PAR(kpp_3d_fields,ncid_landsea,'max_depth',1,1,kpp_3d_fields%ocdepth)
     
     WRITE(6,*) 'Read landsea mask'     
     status=NF_CLOSE(ncid_landsea)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ELSE
     DO ipt=1,npts
        kpp_3d_fields%L_OCEAN(ipt)=.FALSE.
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE mckpp_initialize_landsea
