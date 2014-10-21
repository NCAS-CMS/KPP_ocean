SUBROUTINE mckpp_initialize_relaxation(kpp_3d_fields,kpp_const_fields)
  
! Re-write logic to allow for relaxing either SST or
! salinity - NPK 24/08/11

  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
#include <constants.com>
#include <ocn_advec.com>
#include <couple.com>
#include <sstclim.com>
#include <relax_3d.com>
  
  INTEGER ix,iy,ipoint
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  REAL sst_in(NX_GLOBE,NY_GLOBE,1)
  
  COMMON /save_sstin/ sst_in
  
  DO iy=1,ny
     IF (L_RELAX_SST .and. relax_sst_in(iy) .EQ. 0.0) THEN
        DO ix=1,nx
           ipoint=(iy-1)*nx+ix
           kpp_3d_fields%relax_sst(ipoint)=0.0
        ENDDO
     ELSE
        DO ix=1,nx
           ipoint=(iy-1)*nx+ix
           kpp_3d_fields%relax_sst(ipoint)=1./(relax_sst_in(iy)*kpp_const_fields%spd)
        ENDDO
     ENDIF
     IF (L_RELAX_SAL .and. relax_sal_in(iy) .EQ. 0.0) THEN
        DO ix=1,nx
           ipoint=(iy-1)*nx+ix
           kpp_3d_fields%relax_sal(ipoint)=0.0
        ENDDO
     ELSE
        DO ix=1,nx
           ipoint=(iy-1)*nx+ix
           kpp_3d_fields%relax_sal(ipoint)=1./(relax_sal_in(iy)*kpp_const_fields%spd)
        ENDDO
     ENDIF
     IF (L_RELAX_OCNT .and. relax_ocnt_in(iy) .EQ. 0.0) THEN
        DO ix=1,nx
           ipoint=(iy-1)*nx+ix
           kpp_3d_fields%relax_ocnT(ipoint)=0.0
        ENDDO
     ELSE
        DO ix=1,nx
           ipoint=(iy-1)*nx+ix
           kpp_3d_fields%relax_ocnT(ipoint)=1./(relax_ocnT_in(iy)*kpp_const_fields%spd)
        ENDDO
     ENDIF
  ENDDO
  CALL MCKPP_PHYSICS_OVERRIDES_SST0(kpp_3d_fields)
  DO iy=1,ny
     DO ix=1,nx
        ipoint=(iy-1)*nx+ix
        kpp_3d_fields%fcorr(ipoint)=0.0
        kpp_3d_fields%scorr(ipoint,:)=0.0
     ENDDO
  ENDDO
  
  write(nuout,*) 'calculated SST0, fcorr and scorr'
  
  RETURN
END SUBROUTINE mckpp_initialize_relaxation
