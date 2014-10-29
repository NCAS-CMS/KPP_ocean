SUBROUTINE mckpp_initialize_fluxes_variables(kpp_3d_fields)
  
! Set up parameters for calculating fluxes and initialize fluxes.
! intermediate values computed every ndtld
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
  ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
  TYPE(kpp_3d_type) :: kpp_3d_fields
  integer i,ipt

! extrapolation parameters for fluxes
! Initialize flux arrays
  kpp_3d_fields%wU=0.0
  kpp_3d_fields%wX=0.0
  kpp_3d_fields%wXNT=0.0
  kpp_3d_fields%sflux=0.0

!  call aset(kpp_3d_fields%wU,NZP1tmax*NVP1,0.0)
!  call aset(kpp_3d_fields%wX,NZP1tmax*NSP1,0.0)
!  call aset(kpp_3d_fields%wU,NZP1tmax*NSCLR,0.0)
!  call aset(kpp_3d_fields%sflux,nsflxs*5*(njdt+1),0.0)
  
  DO ipt=1,npts
     do i=1,nsflxs
        kpp_3d_fields%sflux(ipt,i,5,0)=1e-20
     enddo
  ENDDO
  return
end SUBROUTINE mckpp_initialize_fluxes_variables

SUBROUTINE mckpp_initialize_fluxes_file(kpp_const_fields)
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
#include <mc-kpp_3d_type.com>
#include <flx_in.com>
  
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER status,index(3)
  
  index(1)=1
  index(2)=1
  index(3)=1
  
  WRITE(6,*) 'MCKPP_INITIALIZE_FLUXES: Opening file ',kpp_const_fields%forcing_file
  status=NF_OPEN(kpp_const_fields%forcing_file,0,kpp_const_fields%flx_ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'time',kpp_const_fields%flx_timein_id)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'taux',kpp_const_fields%flx_varin_id(1))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'tauy',kpp_const_fields%flx_varin_id(2))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'swf',kpp_const_fields%flx_varin_id(3))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'lwf',kpp_const_fields%flx_varin_id(4))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'lhf',kpp_const_fields%flx_varin_id(5))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'shf',kpp_const_fields%flx_varin_id(6))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'precip',kpp_const_fields%flx_varin_id(7))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  status=NF_GET_VAR1_REAL(kpp_const_fields%flx_ncid,kpp_const_fields%flx_timein_id,&
       index,kpp_const_fields%flx_first_timein)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
            
  RETURN
END SUBROUTINE mckpp_initialize_fluxes_file