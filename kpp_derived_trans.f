      SUBROUTINE kpp_fields_3dto2d(kpp_fields_3d,point,kpp_fields_2d)

      IMPLICIT NONE

c     Accepts a 3D variable of the KPP derived type.
c     Returns a 2D variable of the KPP derived type, extraced from the 3D variable
c     at a specified point.

#include <kpp_3d_type.com>
      TYPE(kpp_3d_type),intent(in)  :: kpp_fields_3d
      TYPE(kpp_2d_type),intent(out) :: kpp_fields_2d
      INTEGER,intent(in) :: point
      INTEGER :: i,j,k
      REAL :: temp
      LOGICAL :: logical_temp

      DO i=1,NZP1         
         DO j=1,NVEL
            temp=kpp_fields_3d%U(point,i,j)
            kpp_fields_2d%U(i,j)=temp
            DO k=0,1
               temp=kpp_fields_3d%Us(point,i,j,k)
               kpp_fields_2d%Us(i,j,k)=temp
            ENDDO
         ENDDO
         DO j=1,NSCLR
            temp=kpp_fields_3d%X(point,i,j)
            kpp_fields_2d%X(i,j)=temp
            DO k=0,1
               temp=kpp_fields_3d%Xs(point,i,j,k)
               kpp_fields_2d%Xs(i,j,k)=temp
            ENDDO
         ENDDO
         temp=kpp_fields_3d%Rig(point,i)
         kpp_fields_2d%Rig(i)=temp
         temp=kpp_fields_3d%Shsq(point,i)
         kpp_fields_2d%Shsq(i)=temp
         temp=kpp_fields_3d%swfrac(point,i)
         kpp_fields_2d%swfrac(i)=temp
         temp=kpp_fields_3d%tinc_fcorr(point,i)
         kpp_fields_2d%tinc_fcorr(i)=temp
         temp=kpp_fields_3d%sinc_fcorr(point,i)
         kpp_fields_2d%sinc_fcorr(i)=temp
         temp=kpp_fields_3d%fcorr_withz(point,i)
         kpp_fields_2d%fcorr_withz(i)=temp      
         temp=kpp_fields_3d%sfcorr_withz(point,i)
         kpp_fields_2d%sfcorr_withz(i)=temp
         temp=kpp_fields_3d%sal_clim(point,i) !Not updated within physics
         kpp_fields_2d%sal_clim(i)=temp
         temp=kpp_fields_3d%ocnT_clim(point,i) !Not updated within physics
         kpp_fields_2d%ocnT_clim(i)=temp
         temp=kpp_fields_3d%ocnTcorr(point,i)
         kpp_fields_2d%ocnTcorr(i)=temp
         IF (i .le. NZ) THEN
            temp=kpp_fields_3d%dbloc(point,i)
            kpp_fields_2d%dbloc(i)=temp
         ENDIF
      ENDDO
      DO i=0,1
         temp=kpp_fields_3d%hmixd(point,i)
         kpp_fields_2d%hmixd(i)=temp
      ENDDO
      DO i=0,NZP1tmax
         temp=kpp_fields_3d%rho(point,i)
         kpp_fields_2d%rho(i)=temp
         temp=kpp_fields_3d%cp(point,i)
         kpp_fields_2d%cp(i)=temp         
         IF (i .gt. 0) THEN
            temp=kpp_fields_3d%buoy(point,i)
            kpp_fields_2d%buoy(i)=temp
         ENDIF
         IF (i .le. NZtmax) THEN
            temp=kpp_fields_3d%difm(point,i)
            kpp_fields_2d%difm(i)=temp
            temp=kpp_fields_3d%difs(point,i)
            kpp_fields_2d%difs(i)=temp
            temp=kpp_fields_3d%dift(point,i)
            kpp_fields_2d%dift(i)=temp
            DO j=1,NVP1
               temp=kpp_fields_3d%wU(point,i,j)
               kpp_fields_2d%wU(i,j)=temp
            ENDDO
            DO j=1,NSP1
               temp=kpp_fields_3d%wX(point,i,j)
               kpp_fields_2d%wX(i,j)=temp
            ENDDO
            DO j=1,NSCLR
               temp=kpp_fields_3d%wXNT(point,i,j)
               kpp_fields_2d%wXNT(i,j)=temp
            ENDDO
            IF (i .gt. 0) THEN
               temp=kpp_fields_3d%ghat(point,i)
               kpp_fields_2d%ghat(i)=temp
            ENDIF
         ENDIF
         IF (i .le. NZ) THEN
            temp=kpp_fields_3d%swdk_opt(point,i)
            kpp_fields_2d%swdk_opt(i)=temp
         ENDIF
      ENDDO
      
      DO i=1,2
         k=kpp_fields_3d%nmodeadv(point,i) ! Not updated within physics
         kpp_fields_2d%nmodeadv(i)=k
         DO j=1,maxmodeadv
            k=kpp_fields_3d%modeadv(point,j,i) !Not updated within physics
            kpp_fields_2d%modeadv(j,i)=k
            temp=kpp_fields_3d%advection(point,j,i) !Not updated within physics
            kpp_fields_2d%advection(j,i)=temp
         ENDDO
      ENDDO

      DO i=1,NSFLXS
         DO j=1,5
            DO k=0,NJDT
               temp=kpp_fields_3d%sflux(point,i,j,k)  !Not updated within physics
               kpp_fields_2d%sflux(i,j,k)=temp
            ENDDO
         ENDDO
      ENDDO

      temp=kpp_fields_3d%ocdepth(point) !Not updated within physics
      kpp_fields_2d%ocdepth=temp      
      logical_temp=kpp_fields_3d%L_OCEAN(point) !Not updated within physics
      kpp_fields_2d%L_OCEAN=logical_temp
      logical_temp=kpp_fields_3d%L_INITFLAG(point)
      kpp_fields_2d%L_INITFLAG=logical_temp
      temp=kpp_fields_3d%f(point) !Not updated within physics
      kpp_fields_2d%f=temp      

      temp=kpp_fields_3d%relax_sst(point) !Not updated within physics
      kpp_fields_2d%relax_sst=temp
      temp=kpp_fields_3d%fcorr(point)
      kpp_fields_2d%fcorr=temp
      temp=kpp_fields_3d%fcorr_twod(point)
      kpp_fields_2d%fcorr_twod=temp
      temp=kpp_fields_3d%SST0(point)      !Not updated within physics
      kpp_fields_2d%SST0=temp
      temp=kpp_fields_3d%relax_sal(point) !Not updated within physics
      kpp_fields_2d%relax_sal=temp
      temp=kpp_fields_3d%relax_ocnt(point) !Not updated within physics
      kpp_fields_2d%relax_ocnt=temp

      temp=kpp_fields_3d%hmix(point)
      kpp_fields_2d%hmix=temp
      temp=kpp_fields_3d%kmix(point)
      kpp_fields_2d%kmix=temp
      temp=kpp_fields_3d%Tref(point)
      kpp_fields_2d%Tref=temp
      temp=kpp_fields_3d%uref(point)
      kpp_fields_2d%uref=temp
      temp=kpp_fields_3d%vref(point)
      kpp_fields_2d%vref=temp
      temp=kpp_fields_3d%Ssurf(point)
      kpp_fields_2d%Ssurf=temp
      temp=kpp_fields_3d%Sref(point) !Not updated within physics
      kpp_fields_2d%Sref=temp
      temp=kpp_fields_3d%SSref(point) !Not updated within physics
      kpp_fields_2d%SSref=temp

      i=kpp_fields_3d%old(point)
      kpp_fields_2d%old=i
      i=kpp_fields_3d%new(point)
      kpp_fields_2d%new=i
      i=kpp_fields_3d%jerlov(point) ! Not updated within physics
      kpp_fields_2d%jerlov=i

      temp=kpp_fields_3d%dlat(point) !Not updated within physics
      kpp_fields_2d%dlat=temp          
      temp=kpp_fields_3d%dlon(point) !Not updated within physics
      kpp_fields_2d%dlon=temp

      RETURN
      END

      SUBROUTINE kpp_fields_2dto3d(kpp_fields_2d,point,kpp_fields_3d)
      IMPLICIT NONE

c     Accepts a 2D and a 3D variable of the KPP derived type.
c     Returns the 3D variable, updated at a specified point with the 
c     values from the 2D variable.

#include <kpp_3d_type.com>
      TYPE(kpp_3d_type),intent(inout)  :: kpp_fields_3d
      TYPE(kpp_2d_type),intent(in) :: kpp_fields_2d
      INTEGER,intent(in) :: point
      INTEGER :: i,j,k
      REAL :: temp
      LOGICAL :: logical_temp

      DO i=1,NZP1
         DO j=1,NVEL
            temp=kpp_fields_2d%U(i,j)
            kpp_fields_3d%U(point,i,j)=temp
            DO k=0,1
               temp=kpp_fields_2d%Us(i,j,k)
               kpp_fields_3d%Us(point,i,j,k)=temp
            ENDDO
         ENDDO
         DO j=1,NSCLR
            temp=kpp_fields_2d%X(i,j)
            kpp_fields_3d%X(point,i,j)=temp
            DO k=0,1
               temp=kpp_fields_2d%Xs(i,j,k)
               kpp_fields_3d%Xs(point,i,j,k)=temp
            ENDDO
         ENDDO
         temp=kpp_fields_2d%Rig(i)
         kpp_fields_3d%Rig(point,i)=temp
         temp=kpp_fields_2d%Shsq(i)
         kpp_fields_3d%Shsq(point,i)=temp
         temp=kpp_fields_2d%swfrac(i)
         kpp_fields_3d%swfrac(point,i)=temp
         temp=kpp_fields_2d%tinc_fcorr(i)
         kpp_fields_3d%tinc_fcorr(point,i)=temp
         temp=kpp_fields_2d%sinc_fcorr(i)
         kpp_fields_3d%sinc_fcorr(point,i)=temp
         temp=kpp_fields_2d%fcorr_withz(i)
         kpp_fields_3d%fcorr_withz(point,i)=temp
         temp=kpp_fields_2d%scorr(i)
         kpp_fields_3d%scorr(point,i)=temp
         temp=kpp_fields_2d%sfcorr_withz(i)
         kpp_fields_3d%sfcorr_withz(point,i)=temp
         temp=kpp_fields_2d%ocnTcorr(i)
         kpp_fields_3d%ocnTcorr(point,i)=temp
         IF (i.le.NZ) THEN 
            temp=kpp_fields_2d%dbloc(i)
            kpp_fields_3d%dbloc(point,i)=temp
         ENDIF
      ENDDO
      DO i=0,1
         temp=kpp_fields_2d%hmixd(i)
         kpp_fields_3d%hmixd(point,i)=temp
      ENDDO
      DO i=0,NZP1tmax
         temp=kpp_fields_2d%rho(i)
         kpp_fields_3d%rho(point,i)=temp
         temp=kpp_fields_2d%cp(i)
         kpp_fields_3d%cp(point,i)=temp         
         IF (i .gt. 0) THEN
            temp=kpp_fields_2d%buoy(i)
            kpp_fields_3d%buoy(point,i)=temp
         ENDIF
         IF (i .le. NZtmax) THEN
            temp=kpp_fields_2d%difm(i)
            kpp_fields_3d%difm(point,i)=temp
            temp=kpp_fields_2d%difs(i)
            kpp_fields_3d%difs(point,i)=temp
            temp=kpp_fields_2d%dift(i)
            kpp_fields_3d%dift(point,i)=temp
            DO j=1,NVP1
               temp=kpp_fields_2d%wU(i,j)
               kpp_fields_3d%wU(point,i,j)=temp
            ENDDO
            DO j=1,NSP1
               temp=kpp_fields_2d%wX(i,j)
               kpp_fields_3d%wX(point,i,j)=temp
            ENDDO
            DO j=1,NSCLR
               temp=kpp_fields_2d%wXNT(i,j)
               kpp_fields_3d%wXNT(point,i,j)=temp
            ENDDO            
            IF (i .gt. 0) THEN 
               temp=kpp_fields_2d%ghat(i)
               kpp_fields_3d%ghat(point,i)=temp
            ENDIF
         ENDIF
         IF (i .le. NZ) THEN
            temp=kpp_fields_2d%swdk_opt(i)
            kpp_fields_3d%swdk_opt(point,i)=temp
         ENDIF
      ENDDO
      
      logical_temp=kpp_fields_2d%L_INITFLAG
      kpp_fields_3d%L_INITFLAG(point)=logical_temp

      temp=kpp_fields_2d%fcorr
      kpp_fields_3d%fcorr(point)=temp
      temp=kpp_fields_2d%fcorr_twod
      kpp_fields_3d%fcorr_twod(point)=temp

      temp=kpp_fields_2d%hmix
      kpp_fields_3d%hmix(point)=temp
      temp=kpp_fields_2d%kmix
      kpp_fields_3d%kmix(point)=temp
      temp=kpp_fields_2d%Tref
      kpp_fields_3d%Tref(point)=temp
      temp=kpp_fields_2d%uref
      kpp_fields_3d%uref(point)=temp
      temp=kpp_fields_2d%vref
      kpp_fields_3d%vref(point)=temp
      temp=kpp_fields_2d%Ssurf
      kpp_fields_3d%Ssurf(point)=temp

      i=kpp_fields_2d%old
      kpp_fields_3d%old(point)=i
      i=kpp_fields_2d%new
      kpp_fields_3d%new(point)=i

      i=kpp_fields_2d%reset_flag
      kpp_fields_3d%reset_flag(point)=i

      RETURN
      END

      SUBROUTINE kpp_const_fields_init(kpp_const_fields)

!     Initialise the const_fields derived type with constants read in from the
!     namelist.  This is a stupid way of doing things, but Fortran won't
!     let me initialise the derived type directly from the namelist.

!     This should be called after *ALL* constants have been read in
!     in steves_3d_ocn.f (subroutine initialize).  It should be called
!     only once.

      IMPLICIT NONE

!     Automatically includes parameter.inc
#include <kpp_3d_type.com>

!     Include all the common blocks containing constants (boo, hiss, common blocks)
#include <constants.com>
#include <flx_paras.com>
#include <ocn_state.com>
#include <ocn_paras.com>
#include <ice_paras.com>
#include <proc_swit.com>
#include <vert_pgrid.com>
#include <timocn.com>
#include <fcorr_in.com>
#include <sfcorr_in.com>
#include <initialcon.com>
#include <ocn_advec.com>
#include <relax_3d.com>

      TYPE(kpp_const_type),intent(inout) :: kpp_const_fields

      kpp_const_fields%spd=spd
      kpp_const_fields%dpy=dpy
      kpp_const_fields%twopi=twopi
      kpp_const_fields%onepi=onepi
      kpp_const_fields%grav=grav
      kpp_const_fields%vonk=vonk
      kpp_const_fields%TK0=TK0
      kpp_const_fields%sbc=sbc
      kpp_const_fields%epsw=epsw
      kpp_const_fields%albocn=albocn
      kpp_const_fields%sice=sice
      kpp_const_fields%EL=EL
      kpp_const_fields%SL=SL
      kpp_const_fields%FL=FL
      kpp_const_fields%FLSN=FLSN

      kpp_const_fields%LKPP=LKPP
      kpp_const_fields%LRI=LRI
      kpp_const_fields%LDD=LDD
      kpp_const_fields%LICE=LICE
      kpp_const_fields%LBIO=LBIO
      kpp_const_fields%LTGRID=LTGRID
      kpp_const_fields%LNBFLX=LNBFLX
      kpp_const_fields%LRHS=LRHS
      kpp_const_fields%L_SSref=L_SSref
      kpp_const_fields%L_RELAX_SST=L_RELAX_SST
      kpp_const_fields%L_RELAX_CALCONLY=L_RELAX_CALCONLY
      kpp_const_fields%L_FCORR=L_FCORR
      kpp_const_fields%L_SFCORR=L_SFCORR
      kpp_const_fields%L_FCORR_WITHZ=L_FCORR_WITHZ
      kpp_const_fields%L_SFCORR_WITHZ=L_SFCORR_WITHZ
      kpp_const_fields%L_RESTART=L_RESTART
      kpp_const_fields%L_RELAX_SAL=L_RELAX_SAL
      kpp_const_fields%L_RELAX_OCNT=L_RELAX_OCNT

      RETURN
      END
      
      
      
