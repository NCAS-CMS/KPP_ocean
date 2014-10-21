SUBROUTINE mckpp_physics_overrides_bottomtemp(kpp_3d_fields,kpp_const_fields)  

  ! Written by NPK 10/4/08
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <mc-kpp_3d_type.com>
  
  INTEGER ipt,z
  REAL bottom_temp(NPTS)
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  DO ipt=1,npts
     kpp_3d_fields%tinc_fcorr(ipt,NZP1)=kpp_3d_fields%bottom_temp(ipt)-kpp_3D_fields%X(ipt,NZP1,1)
     kpp_3d_fields%ocnTcorr(ipt,NZP1)=kpp_3d_fields%tinc_fcorr(ipt,NZP1)*&
          kpp_3d_fields%rho(ipt,NZP1)*kpp_3d_fields%cp(ipt,NZP1)/&
          kpp_const_fields%dto
     kpp_3d_fields%X(ipt,NZP1,1) = kpp_3d_fields%bottom_temp(ipt)
  ENDDO
  
  RETURN
END SUBROUTINE mckpp_physics_overrides_bottomtemp

SUBROUTINE mckpp_physics_overrides_sst0(kpp_3d_fields)
  ! Written by NPK 27/8/07
  
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
  !include 'constants.com'
  !include 'ocn_advec.com'
#include <couple.com>
#include <sstclim.com>

  INTEGER ix,iy,ipoint
  
  TYPE(kpp_3d_type) :: kpp_3d_fields
  REAL sst_in(NX_GLOBE,NY_GLOBE,1)
  COMMON /save_sstin/ sst_in
  
  DO iy=1,ny
     DO ix=1,nx
        ipoint=(iy-1)*nx+ix
        kpp_3d_fields%SST0(ipoint)=SST_in(ix+ifirst-1,iy+jfirst-1,1)
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE mckpp_physics_overrides_sst0

SUBROUTINE mckpp_physics_overrides_check_profile(kpp_1d_fields,kpp_const_fields)
  IMPLICIT NONE

#include <mc-kpp_3d_type.com>

  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  INTEGER :: z,j
  REAL :: dz_total,dtdz_total,dz

  ! If the integration has failed because of unrealistic values in T, S, U or V
  ! or very high RMS difference between the old and new profiles, then reset
  ! T and S to climatology (if available) and U and V to the initial profiles.
  ! NPK 17/5/13.
  IF (kpp_1d_fields%comp_flag .and. kpp_const_fields%ocnT_file .ne. 'none' .and. &
       kpp_const_fields%sal_file .ne. 'none') THEN
     WRITE(6,*) 'Resetting point to climatology ...'
     kpp_1d_fields%X(:,1)=kpp_1d_fields%ocnT_clim(:)
     WRITE(6,*) 'T = ',kpp_1d_fields%ocnT_clim(:)
     kpp_1d_fields%X(:,2)=kpp_1d_fields%sal_clim(:)
     WRITE(6,*) 'S = ',kpp_1d_fields%sal_clim(:)
     kpp_1d_fields%U=kpp_1d_fields%U_init(:,:)
     WRITE(6,*) 'U = ',kpp_1d_fields%U_init(:,1)
     WRITE(6,*) 'V = ',kpp_1d_fields%U_init(:,2)
     kpp_1d_fields%reset_flag=999
  ELSE IF (kpp_1d_fields%comp_flag) THEN
     WRITE(6,*) 'Cannot reset point to T,S climatology as either ocean temperature or salinity data '//&
          'not provided.  Will reset currents to initial conditions and keep going.'
     kpp_1d_fields%U=kpp_1d_fields%U_init(:,:)
     kpp_1d_fields%reset_flag=999
  ENDIF
  
  ! Check whether the temperature at any (x,z) point is less than the
  ! threshold for sea ice (-1.8C).  If it is, reset it to -1.8C and
  ! set a flag.  The flag can be requested as a diagnostic (singout 9).
  ! Note that the value of the flag is equal to the *fraction* of levels
  ! at that point that were < -1.8C.
  IF (kpp_1d_fields%L_OCEAN) THEN
     DO z=1,NZP1
        IF (kpp_1d_fields%X(z,1) .lt. -1.8) THEN
           kpp_1d_fields%tinc_fcorr(z)=kpp_1d_fields%tinc_fcorr(z)+&
                (-1.8-kpp_1d_fields%X(z,1))
           kpp_1d_fields%X(z,1)=-1.8
           kpp_1d_fields%freeze_flag=kpp_1d_fields%freeze_flag+1.0/REAL(NZP1)
        ENDIF
     ENDDO
  ENDIF
        
  ! Check whether the temperature difference between the surface
  ! and a user-specified level (presumably deep) is less than a user-specified 
  ! threshold (presumably small).  If so, reset the temperature and salinity
  ! profiles to climatological values.  Added to prevent spurious very
  ! deep mixing that creates unrealistic isothermal (and isohaline) layers.
  ! NPK 15/5/2013 for R4.
  IF (kpp_1d_fields%L_OCEAN) THEN
     dtdz_total=0.
     dz_total=0.
     DO j=2,kpp_const_fields%iso_bot
        dz=kpp_const_fields%zm(j)-kpp_const_fields%zm(j-1)
        dtdz_total=dtdz_total+ABS((kpp_1d_fields%X(j,1)-&
             kpp_1d_fields%X(j-1,1)))*dz
        dz_total=dz_total+dz
     ENDDO
     dtdz_total=dtdz_total/dz_total
     
     ! If resetting to climatology because of isothermal layer (rather than because of 
     ! computational instability trap in ocn.f), then set reset_flag to a negative
     ! value (-1*number of interations in of semi-implicit integration in ocn.f).
     IF (ABS(dtdz_total).lt.kpp_const_fields%iso_thresh) THEN
        kpp_1d_fields%X(:,1)=kpp_1d_fields%ocnT_clim(:)
        kpp_1d_fields%X(:,2)=kpp_1d_fields%sal_clim(:)
        kpp_1d_fields%reset_flag=(-1.)*kpp_1d_fields%reset_flag
     ENDIF
  ELSE
     kpp_1d_fields%reset_flag=0
  ENDIF
  
  RETURN
END SUBROUTINE mckpp_physics_overrides_check_profile

