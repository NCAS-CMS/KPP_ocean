      PROGRAM ocn_model_3D
**************************************************************************
* 3D version of the 1D ocean model using the kpp mixing scheme of
* Large et al, with his interface between the ocean and kpp scheme
* calls his subroutines
* init_ocn : to initialize the ocean model
* ocn_step : to update the model
*
* Also uses his parameter namelists and common blocks
* There maybe a lot of unnecessary variables and parameters here, but until
* we get something working it seems foolish to try and strip too much out.
* Written April 2002
*
*  Steve Woolnough
**************************************************************************

c      USE kpp_type_mod
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
#include <landsea.com>

#ifdef COUPLE
#ifdef OASIS2
#include <param.h>
#endif /*OASIS2*/
#ifdef OASIS3
#include <kpp_oasis3.inc>
#endif /*OASIS3*/
#endif /*COUPLE*/

#include <bottomclim.com>
#include <currclim.com>
#include <couple.com>
#include <fcorr_in.com>
#include <sfcorr_in.com>
#include <relax_3d.com>
#include <times.com>
#include <timocn.com>
#include <vert_pgrid.com>
#include <local_pt.com>
#include <output.com>
#include <initialcon.com>
#include <sstclim.com>
#include <ocn_advec.com>

* Local variables
      TYPE(kpp_3d_type),allocatable :: kpp_3d_fields
      TYPE(kpp_2d_type) :: kpp_2d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      TYPE(kpp_timer_type) :: kpp_timer
      REAL,allocatable :: VEC_mean(:,:,:),SCLR_mean(:,:),bottom_temp(:),
     +     VEC_range(:,:,:,:),SCLR_range(:,:,:)

      INTEGER j,k,nflx,nstep,ix,iy,ipt_globe

c Initialize the 3D KPP model
c Setup the constants, read the namelists, setup the initial conditions
c
c Add bottom temperatures to the call to INITIALIZE, in case the bottom temperatures
c are initialized but never updated.  Previously, they never would have made it back
c into the main program.  NPK 17/08/10 - R3

      INTEGER (KIND=2) :: nthreads,tid,omp_get_num_threads,
     +  omp_get_thread_num
      CHARACTER(LEN=21) phys_timer_name
      CHARACTER(LEN=19) trans_timer_name

      CALL KPP_TIMER_INIT(kpp_timer)
#ifdef OPENMP
!$OMP PARALLEL PRIVATE(nthreads)
      CALL OMP_SET_DYNAMIC(.FALSE.)
      nthreads=OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
      WRITE(6,*) 'Initialising ',nthreads,'timers'
      IF (nthreads .eq. 0 .or. nthreads .ge. 100) THEN
         WRITE(6,*) 'nthreads = 0, resetting nthreads = 24'
         nthreads=24
         CALL OMP_SET_NUM_THREADS(24)
      ENDIF
      DO k=0,nthreads-1
        WRITE(phys_timer_name,'(A19,I2)') 'KPP Physics thread ',k
        CALL KPP_TIMER_TIME(kpp_timer,phys_timer_name,1)
        CALL KPP_TIMER_TIME(kpp_timer,phys_timer_name,0)
        WRITE(trans_timer_name,'(A17,I2)') 'KPP 3D/2D thread ',k
        CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,1)
        CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,0)
      ENDDO
#endif
      allocate(VEC_mean(NPTS,NZP1,NVEC_MEAN))
      allocate(SCLR_mean(NPTS,NSCLR_MEAN))
      allocate(VEC_range(NPTS,NZP1,NVEC_RANGE,2))
      allocate(SCLR_range(NPTS,NSCLR_RANGE,2))
      allocate(bottom_temp(NPTS))
      allocate(kpp_const_fields%wmt(0:891,0:49))
      allocate(kpp_const_fields%wst(0:891,0:49))
      allocate(kpp_const_fields%tri(0:NZtmax,0:1,NGRID))
      allocate(kpp_3d_fields)

      CALL KPP_TIMER_TIME(kpp_timer,'Initialize',1)
      CALL initialize(kpp_3d_fields,kpp_const_fields,bottom_temp,
     +     VEC_mean,SCLR_mean,VEC_range,SCLR_range)
      CALL KPP_TIMER_TIME(kpp_timer,'Initialize',0)
      CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)

      write(nuout,*) '3D KPP model initialized with, ',NZ,' levels'
      write(nuout,*) 'and ',npts,' points'

#ifdef COUPLE
#ifdef OASIS2
      IF ((im .NE. NX_GLOBE) .OR. (jm .NE. NY_GLOBE)) THEN
c
c     Test that KPP and OASIS have the same horizontal grid sizes.
c     These are controlled by NX_GLOBE and NY_GLOBE in parameter.inc
c     (for KPP) and by im and jm in param.h (for OASIS).
c     This applies to OASIS2 only.
c
         WRITE(nuout,*) 'KPP : KPP and OASIS2 do not have the same ',
     +        'horizontal grid sizes.  KPP has NX_GLOBE=',NX_GLOBE,
     +        'and NY_GLOBE=',NY_GLOBE,'while OASIS2 has im=',im,
     +        'and jm=',jm,'.  You can control the KPP settings in ',
     +        'parameter.inc and the OASIS2 settings in param.h.'
         CALL halte('im,jm not equal to NX_GLOBE,NY_GLOBE')
      ELSE
c     Initialize the OASIS2 coupling interface
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
         CALL KPP_TIMER_TIME(kpp_timer,'OASIS2 initialization',1)
         CALL inicmo(nend*ndtocn,ndtocn,int(kpp_const_fields%dto))
         CALL KPP_TIMER_TIME(kpp_timer,'OASIS2 initialization',0)
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
      ENDIF
#endif /*OASIS2*/
#ifdef OASIS3
      CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
      CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 initialization',1)
      CALL mpi1_oasis3_init()
      CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 initialization',0)
      CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
c
c     When coupling to the UM via OASIS3, we must send the initial ocean fields to
c     the atmosphere.  The UM does an OASIS3 "get" first, then an OASIS3 "send."
c     We must match this with a "send" before we post our first "get", otherwise
c     the whole system deadlocks and no one has any fun.
c
c     Note that the OASIS3 "get" is handled by the first call to <fluxes> in
c     the time-stepping loop below.     NPK 2/10/09
c
      CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
      CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 output',1)
      CALL mpi1_oasis3_output(kpp_3d_fields,kpp_const_fields)
      CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 output',0)
      CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
#endif /*OASIS3*/
#endif /*COUPLE*/

      WRITE(nuout,*) ' Starting the ocean integration'
c
c     Main KPP time-stepping loop
c
      DO ntime=1,nend*ndtocn
         kpp_const_fields%ntime=ntime
         IF (MOD(kpp_const_fields%ntime-1,ndtocn) .EQ. 0) THEN
!            WRITE(6,*) 'Calling fluxes'
#ifdef COUPLE
#ifdef OASIS2
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'OASIS2 input and MPI wait',1)
            CALL fluxes(kpp_3d_fields,kpp_const_fields,kpp_timer)
            CALL KPP_TIMER_TIME(kpp_timer,'OASIS2 input and MPI wait',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
#endif /*OASIS2*/
#ifdef OASIS3
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 input and MPI wait',1)
            CALL fluxes(kpp_3d_fields,kpp_const_fields,kpp_timer)
            CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 input and MPI wait',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
#endif /*OASIS3*/
#else
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Update surface fluxes',1)
            CALL fluxes(kpp_3d_fields,kpp_const_fields,kpp_timer)
            CALL KPP_TIMER_TIME(kpp_timer,'Update surface fluxes',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
#endif /*COUPLE*/
!            WRITE(6,*) 'Called fluxes'
         ENDIF
c
c     Re-writing update logic to allow user to choose independently whether to
c     update SST, sea ice, flux corrections and bottom temperatures.
c     NPK 15/10/09 - R3
c
         IF (L_UPD_CLIMSST .AND. MOD(ntime-1,ndtupdsst) .EQ. 0) THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
            CALL read_sstin(kpp_3d_fields,kpp_const_fields)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
            WRITE(nuout,*) 'KPP: Called read_sstin, ntime =',
     +           kpp_const_fields%ntime
c
c     SST0 specifies the temperature to which to relax the mixed-layer temperature.
c     The use of the relaxation is controlled by L_RELAX_SST.
c     Changing L_RELAX_SST .OR. L_COUPLE to just L_RELAX_SST - NPK 2/11/09 - R3
c
            IF (L_RELAX_SST) THEN
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
               CALL KPP_TIMER_TIME(kpp_timer,'SST relaxation',1)
               CALL upd_sst0(kpp_3d_fields)
               CALL KPP_TIMER_TIME(kpp_timer,'SST relaxation',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
               WRITE(nuout,*) 'KPP: Called upd_sst0, ntime =',
     +              kpp_const_fields%ntime
            ENDIF
         ENDIF
         IF (L_UPD_CLIMICE .AND. MOD(ntime-1,ndtupdice) .EQ. 0) THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
            CALL read_icein(kpp_3d_fields,kpp_const_fields)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
            WRITE(nuout,*) 'KPP: Called read_icein, ntime =',
     +           kpp_const_fields%ntime
         ENDIF
         IF (L_UPD_CLIMCURR .AND. MOD(ntime-1,ndtupdcurr) .EQ. 0) THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
            CALL read_surface_currents(kpp_3d_fields,kpp_const_fields)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
            WRITE(nuout,*) 'KPP: Called read_surface_currents, ntime =',
     +           kpp_const_fields%ntime
         ENDIF
         IF (L_UPD_FCORR .AND. MOD(ntime-1,ndtupdfcorr) .EQ. 0) THEN
            IF (L_FCORR_WITHZ) THEN
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
               CALL read_fcorrwithz(kpp_3d_fields,kpp_const_fields)
               CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
               WRITE(nuout,*) 'KPP: Called read_fcorrwithz, ntime =',
     +              kpp_const_fields%ntime
            ELSEIF (L_FCORR) THEN
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
               CALL read_fcorr(kpp_3d_fields,kpp_const_fields)
               CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
               WRITE(nuout,*) 'KPP: Called read_fcorr, ntime =',
     +              kpp_const_fields%ntime
            ENDIF
         ENDIF
!added SFCORR LH 24/05/2013
         IF (L_UPD_SFCORR .AND. MOD(ntime-1,ndtupdsfcorr) .EQ. 0) THEN
            IF (L_SFCORR_WITHZ) THEN
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
               CALL read_sfcorrwithz(kpp_3d_fields,kpp_const_fields)
               CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
               WRITE(nuout,*) 'KPP: Called read_sfcorrwithz, ntime =',
     +              kpp_const_fields%ntime
            ELSEIF (L_SFCORR) THEN
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
               CALL read_sfcorr(kpp_3d_fields,kpp_const_fields)
               CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
               WRITE(nuout,*) 'KPP: Called read_sfcorr, ntime =',
     +              kpp_const_fields%ntime
            ENDIF
         ENDIF
         IF (L_UPD_BOTTOM_TEMP .AND. MOD(ntime-1,ndtupdbottom) .EQ. 0)
     +        THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
            CALL read_bottom_temp(kpp_3d_fields,kpp_const_fields,
     +           bottom_temp)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
            WRITE(nuout,*) 'KPP: Called read_bottom_temp, ntime =',
     &           kpp_const_fields%ntime
         ENDIF
c
c     Add ability to relax to salinity climatology
c     NPK 24/08/11 - R3b3
         IF (L_UPD_SAL .AND. MOD(ntime-1,ndtupdsal) .EQ. 0 .AND.
     +        .NOT. L_INTERP_SAL) THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
            CALL read_salinity(kpp_3d_fields,kpp_const_fields)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
            WRITE(nuout,*) 'KPP: Called read_salinity, ntime =',
     +           kpp_const_fields%ntime
         ELSEIF (L_UPD_SAL .AND. L_INTERP_SAL .AND.
     +           MOD(ntime-1,ndt_interp_sal).EQ.0) THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
            CALL INTERP_SAL(kpp_3d_fields,kpp_const_fields)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
            WRITE(nuout,*) 'KPP: Interpolated ocean salinity,'//
     +           ' ntime =',kpp_const_fields%ntime
         ENDIF
         IF (L_UPD_OCNT .AND. MOD(ntime-1,ndtupdocnt) .EQ. 0
     +        .AND. .NOT. L_INTERP_OCNT) THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
            CALL read_ocean_temperatures(kpp_3d_fields,kpp_const_fields)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
            WRITE(nuout,*) 'KPP: Called read_ocean_temperatures,'//
     +           ' ntime =',kpp_const_fields%ntime
         ELSEIF (L_UPD_OCNT .AND. L_INTERP_OCNT .AND.
     +           MOD(ntime-1,ndt_interp_ocnt).EQ.0) THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
            CALL INTERP_OCNT(kpp_3d_fields,kpp_const_fields)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
            WRITE(nuout,*) 'KPP: Interpolated ocean temperatures,'//
     +           ' ntime =',kpp_const_fields%ntime
         ENDIF
c
c     Call Large et al. (1994) boundary-layer ocean physics routines
c
         kpp_const_fields%time=kpp_const_fields%startt+ntime*
     +        kpp_const_fields%dto/kpp_const_fields%spd
!         WRITE(nuout,*) 'KPP: Entering ocnstep loop, ntime=',ntime
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
c         CALL KPP_TIMER_TIME(kpp_timer,'KPP Physics (all)',1)
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(kpp_3d_fields,kpp_const_fields)
!$OMP& SHARED(kpp_timer,ocnT_file,sal_file,ifirst,jfirst)
!$OMP& PRIVATE(trans_timer_name,phys_timer_name,tid)
         tid=OMP_GET_THREAD_NUM()
         WRITE(trans_timer_name,'(A17,I2)') 'KPP 3D/2D thread ',tid
         WRITE(phys_timer_name,'(A19,I2)') 'KPP Physics thread ',tid
!$OMP DO SCHEDULE(dynamic)
#else
         WRITE(trans_timer_name,'(A19)') 'KPP 3D/2D thread 01'
         WRITE(phys_timer_name,'(A21)') 'KPP Physics thread 01'
#endif
#ifdef COUPLE
         DO ix=1,NX_GLOBE
            DO iy=1,NY_GLOBE
               ipt_globe=(iy-1)*NX_GLOBE+ix
               ipt=(iy-jfirst)*NX+(ix-ifirst)+1
               IF (kpp_3d_fields%L_OCEAN(ipt) .and.
     +              kpp_3d_fields%cplwght(ipt_globe) .gt. 0) THEN
#else
        DO ipt=1,npts
!           IF (iy .ge. jfirst .and. iy .le. jlast .and.
!     +                 ix .ge. ifirst .and. ix .le. ilast) THEN
           !ipt=(iy-jfirst)*NX+(ix-ifirst)+1
           IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
!                        WRITE(6,*) 'ipt = ',ipt,'ix=',ix,'iy=',iy
#endif
                  CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,1)
                  CALL kpp_fields_3dto2d(kpp_3d_fields,ipt,
     +                 kpp_2d_fields)
                  CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,0)
                  CALL KPP_TIMER_TIME(kpp_timer,phys_timer_name,1)
                  CALL ocnstep(kpp_2d_fields,kpp_const_fields)
c     If the integration has failed because of unrealistic values in T, S, U or V
c     or very high RMS difference between the old and new profiles, then reset
c     T and S to climatology (if available) and U and V to the initial profiles.
c     NPK 17/5/13.
                  IF (kpp_2d_fields%comp_flag .and. ocnT_file .ne.
     +                 'none' .and. sal_file .ne. 'none') THEN
                     WRITE(6,*) 'Resetting point to climatology ...'
                     kpp_2d_fields%X(:,1)=kpp_2d_fields%ocnT_clim(:)
                     !WRITE(6,*) 'T = ',kpp_2d_fields%ocnT_clim(:)
                     kpp_2d_fields%X(:,2)=kpp_2d_fields%sal_clim(:)
                     !WRITE(6,*) 'S = ',kpp_2d_fields%sal_clim(:)
                     kpp_2d_fields%U=kpp_3d_fields%U_init(ipt,:,:)
                     !WRITE(6,*) 'U = ',kpp_3d_fields%U_init(ipt,:,1)
                     !WRITE(6,*) 'V = ',kpp_3d_fields%U_init(ipt,:,2)
                     kpp_2d_fields%reset_flag=999
                  ELSE IF (kpp_2d_fields%comp_flag) THEN
                     WRITE(6,*) 'Cannot reset pt to T,S climatology '//
     +                    'as either ocean T or S data '//
     +                    'not provided.  Will reset U,V to initial'//
     +                    'conditions and keep going.'
                     kpp_2d_fields%U=kpp_3d_fields%U_init(ipt,:,:)
                     kpp_2d_fields%reset_flag=999
                  ENDIF
                  CALL KPP_TIMER_TIME(kpp_timer,phys_timer_name,0)
                  CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,1)
                  CALL kpp_fields_2dto3d(kpp_2d_fields,ipt,
     +                 kpp_3d_fields)
                  CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,0)
               ENDIF
!#ifndef COUPLE
!            ENDIF
!#endif
#ifdef COUPLE
            ENDDO
#endif
         ENDDO
#ifdef OPENMP
!$OMP END DO
c	 WRITE(6,*) 'OpenMP thread ',tid,' finished main DO loop'
!$OMP END PARALLEL
#endif
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
c
c     Following the physics routines, update the temperature of the
c     bottom layer, if necessary
c     NPK 10/4/08 - R1
c
         IF (L_VARY_BOTTOM_TEMP) THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
            CALL upd_bottom_temp(kpp_3d_fields,kpp_const_fields,
     +           bottom_temp)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
         ENDIF
         IF (L_NO_FREEZE) THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
            kpp_3d_fields%freeze_flag(:)=0.
            CALL check_freezing(kpp_3d_fields)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
         ENDIF
         IF (L_NO_ISOTHERM) THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
            CALL check_isothermal(kpp_3d_fields,kpp_const_fields)
            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
         ENDIF
c         IF (L_DAMP_CURR) THEN
c            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
c            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
c            kpp_3d_fields%dampu_flag(:)=0.
c            kpp_3d_fields%dampv_flag(:)=0.
c            CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
c            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
c         ENDIF

!     If using a slab ocean, enforce constant salinity to prevent large
!     drifts in long simulations that can affect heat capacity and density.
!     Could make this controllable from a namelist option, but precise value
!     used probably doesn't matter too much - just need to prevent highly unrealistic
!     values developing due to lack of corrections.
!     NPK 13/09/16
         IF (kpp_const_fields%L_SLAB) THEN
!     Values are held as deviations from reference salinity, so set prognostic to zero
!     to preserve initial conditions.
            kpp_3d_fields%X(:,1,2) = 0.0
         ENDIF            

!         WRITE(nuout,*) 'KPP: Finished ocnstep loop, ntime=',ntime
c
c     Implement more-frequent checkpointing, upon request
c     NPK 02/02/10 - R3
c
         IF (ndt_per_restart .NE. nend*ndtocn) THEN
            IF (MOD(ntime,ndt_per_restart).EQ.0) THEN
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Writing Restart File',1)
               IF (kpp_const_fields%time .lt. 10) THEN
                  WRITE(restart_time,'(A4,I1)') '0000',
     +                 FLOOR(kpp_const_fields%time)
               ELSEIF (kpp_const_fields%time .lt. 100) THEN
                  WRITE(restart_time,'(A3,I2)') '000',
     +                 FLOOR(kpp_const_fields%time)
               ELSEIF (kpp_const_fields%time .lt. 1000) THEN
                  WRITE(restart_time,'(A2,I3)') '00',
     +                 FLOOR(kpp_const_fields%time)
               ELSEIF (kpp_const_fields%time .lt. 10000) THEN
                  WRITE(restart_time,'(A1,I4)') '0',
     +                 FLOOR(kpp_const_fields%time)
               ELSE
                  WRITE(restart_time,'(I5)')
     +                 FLOOR(kpp_const_fields%time)
               ENDIF
               WRITE(restart_outfile,'(A12,A5)') 'KPP.restart.',
     +              restart_time
               CALL WRITE_RESTART(kpp_3d_fields,kpp_const_fields,
     +              restart_outfile)
               CALL KPP_TIMER_TIME(kpp_timer,'Writing Restart File',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
            ENDIF
         ENDIF
!         WRITE(nuout,*) 'KPP: Finished checkpointing'
c
c     If we are going to output means, take means at the end of the timestep
c     NPK 25/2/08 - R1
c
         IF (L_OUTPUT_MEAN) THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Time-meaning of output',1)
            CALL mean_output(kpp_3d_fields,VEC_mean,SCLR_mean)
!            WRITE(6,*) 'KPP: Finished meaning output'
            CALL KPP_TIMER_TIME(kpp_timer,'Time-meaning of output',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
         ENDIF
!         WRITE(6,*) 'KPP: Finished meaning output'
         IF (L_OUTPUT_RANGE) THEN
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Time-meaning of output',1)
            CALL range_output(kpp_3d_fields,VEC_range,SCLR_range)
            CALL KPP_TIMER_TIME(kpp_timer,'Time-meaning of output',0)
            CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
         ENDIF
!         WRITE(nuout,*) 'KPP: Finished ranging output'
c
c     Output means every ndtout_mean timesteps
c     NPK 25/2/08 - R1
c
         DO j=1,N_VAROUTS
            IF (ndt_varout_mean(j) .gt. 0 .and. L_OUTPUT_MEAN) THEN
               IF (MOD(ntime,ndt_varout_mean(j)).EQ.0) THEN
                  CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
                  CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
                  CALL write_means(kpp_3d_fields,kpp_const_fields,
     +                 VEC_mean,SCLR_mean,j,varid_vec_mean(j),
     +                 zprof_varout_mean(j),ntout_vec_mean(j))
                  CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
                  CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
               ENDIF
            ENDIF
            IF (ndt_varout_inst(j) .gt. 0 .and. L_OUTPUT_INST) THEN
               IF (MOD(ntime,ndt_varout_inst(j)).EQ.0) THEN
                  CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
                  CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
                  CALL output_inst(kpp_3d_fields,kpp_const_fields,j,
     +                 varid_vec(j),zprof_varout_inst(j),
     +                 ntout_vec_inst(j))
                  CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
                  CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
               ENDIF
            ENDIF
            IF (ndt_varout_range(j) .gt. 0 .and. L_OUTPUT_RANGE) THEN
               IF (MOD(ntime,ndt_varout_range(j)).EQ.0) THEN
                  CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
                  CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
                  CALL output_range(kpp_3d_fields,kpp_const_fields,
     +                 VEC_range,SCLR_range,j,
     +                 varid_vec_range(j),zprof_varout_range(j),
     +                 ntout_vec_range(j))
                  CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
                  CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
               ENDIF
            ENDIF
         ENDDO
         DO j=1,N_SINGOUTS
            IF (ndt_singout_mean(j) .gt. 0 .and. L_OUTPUT_MEAN) THEN
               IF (MOD(ntime,ndt_singout_mean(j)).EQ.0) THEN
                  CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
                  CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
                  CALL write_means(kpp_3d_fields,kpp_const_fields,
     +                 VEC_mean,SCLR_mean,j+N_VAROUTS,varid_sing_mean(j)
     +                 ,0,ntout_sing_mean(j))
                  CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
                  CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
               ENDIF
            ENDIF
            IF (ndt_singout_inst(j) .gt. 0 .and. L_OUTPUT_INST) THEN
               IF (MOD(ntime,ndt_singout_inst(j)).EQ.0) THEN
                  CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
                  CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
                  CALL output_inst(kpp_3d_fields,kpp_const_fields,
     +                 j+N_VAROUTS,varid_sing(j),0,ntout_sing_inst(j))
                  CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
                  CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
               ENDIF
            ENDIF
            IF (ndt_singout_range(j) .gt. 0 .and. L_OUTPUT_RANGE) THEN
               IF (MOD(ntime,ndt_singout_range(j)) .EQ. 0) THEN
                  CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
                  CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
                  CALL output_range(kpp_3d_fields,kpp_const_fields,
     +                 VEC_range,SCLR_range,j+N_VAROUTS,
     +                 varid_sing_range(j),0,ntout_sing_range(j))
                  CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
                  CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
               ENDIF
            ENDIF
         ENDDO
c
c     If we need to make new instantaneous and mean output files,
c     then close the current files and create new ones.
c     Updated to take advantage of new ndt_per_file option.
c     NPK 10/6/09 - R2
c
         IF (kpp_const_fields%time .GE. float(day_out) .AND. ntime
     +        .NE. nend*ndtocn) THEN
            day_out=day_out+NINT(kpp_const_fields%dtsec/FLOAT(ndtocn)
     +           *FLOAT(ndt_per_file)/kpp_const_fields%spd)
            !WRITE(nuout,*) 'day_out=',day_out
            IF (L_OUTPUT_INST) THEN
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
               write(output_file(flen+2:flen+6),'(i5.5)') day_out
               CALL output_close(ncid_out)
	       ntout_vec_inst(:)=1
	       ntout_sing_inst(:)=1
               CALL init_output(output_file,ncid_out,
     +              kpp_3d_fields,kpp_const_fields,ndt_varout_inst,
     +              ndt_singout_inst,varid_vec,varid_sing,
     +              zprof_varout_inst,.FALSE.,.TRUE.)
               CALL output_open(output_file,ncid_out)
               CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
            ENDIF
            IF (L_OUTPUT_MEAN) THEN
               write(mean_output_file(flen+2:flen+6),'(i5.5)')
     &              day_out
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
               CALL output_close(mean_ncid_out)
	       ntout_vec_mean(:)=1
	       ntout_sing_mean(:)=1
               CALL init_output(mean_output_file,mean_ncid_out,
     +              kpp_3d_fields,kpp_const_fields,ndt_varout_mean,
     +              ndt_singout_mean,varid_vec_mean,varid_sing_mean,
     +              zprof_varout_mean,.TRUE.,.FALSE.)
               CALL output_open(mean_output_file,mean_ncid_out)
               nout_mean=1
               CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
            ENDIF
            IF (L_OUTPUT_RANGE) THEN
               write(min_output_file(flen+2:flen+6),'(i5.5)')
     &              day_out
               write(max_output_file(flen+2:flen+6),'(i5.5)')
     &              day_out
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
               CALL output_close(min_ncid_out)
               CALL output_close(max_ncid_out)
	       ntout_vec_range(:)=1
	       ntout_sing_range(:)=1
               CALL init_output(min_output_file,min_ncid_out,
     +              kpp_3d_fields,kpp_const_fields,ndt_varout_range,
     +              ndt_singout_range,varid_vec_range,varid_sing_range,
     +              zprof_varout_range,.FALSE.,.FALSE.)
               CALL init_output(max_output_file,max_ncid_out,
     +              kpp_3d_fields,kpp_const_fields,ndt_varout_range,
     +              ndt_singout_range,varid_vec_range,varid_sing_range,
     +              zprof_varout_range,.FALSE.,.FALSE.)
               CALL output_open(min_output_file,min_ncid_out)
               CALL output_open(max_output_file,max_ncid_out)
               CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
               CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
            ENDIF
         ENDIF

c
c     Output coupled fields to the atmospheric models via
c     the approprite coupling routines.
c     NPK May 2008 for OASIS2 - R1
c     NPK March 2009 for CFS - R2
c     NPK 18/09/09 for OASIS3 - R3
c
#ifdef COUPLE
         IF ((MOD(ntime,ndtocn) .EQ. 0)
#ifdef OASIS2
     +        .AND. ntime .NE. nend*ndtocn) THEN
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
         CALL KPP_TIMER_TIME(kpp_timer,'OASIS2 output',1)
         CALL coupled_out(kpp_3d_fields%X,ntime,.FALSE.)
         CALL KPP_TIMER_TIME(kpp_timer,'OASIS2 output',0)
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
#else
#ifdef OASIS3
     +        .AND. ntime .NE. nend*ndtocn) THEN
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
         CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 output',1)
         CALL mpi1_oasis3_output(kpp_3d_fields,kpp_const_fields)
         CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 output',0)
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
#else
#ifdef CFS
     +        ) THEN
            CALL CFS_WRITE_GRIB_SSTS(kpp_3d_fields)
#endif /*CFS*/
#endif /*OASIS3*/
#endif /*OASIS2*/
         ENDIF
#endif /*COUPLE*/
      ENDDO

      WRITE(nuout,*) 'KPP: Successful termination of the model ',
     +     'integration'

      IF (L_RESTARTW) THEN
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
         CALL KPP_TIMER_TIME(kpp_timer,'Writing restart file',1)
         IF (kpp_const_fields%time .lt. 10) THEN
            WRITE(restart_time,'(A4,I1)') '0000',
     +           FLOOR(kpp_const_fields%time)
         ELSEIF (kpp_const_fields%time .lt. 100) THEN
            WRITE(restart_time,'(A3,I2)') '000',
     +           FLOOR(kpp_const_fields%time)
         ELSEIF (kpp_const_fields%time .lt. 1000) THEN
            WRITE(restart_time,'(A2,I3)') '00',
     +           FLOOR(kpp_const_fields%time)
         ELSEIF (kpp_const_fields%time .lt. 10000) THEN
            WRITE(restart_time,'(A1,I4)') '0',
     +           FLOOR(kpp_const_fields%time)
         ELSE
            WRITE(restart_time,'(I5)')
     +           FLOOR(kpp_const_fields%time)
         ENDIF
         WRITE(restart_outfile,'(A12,A5)') 'KPP.restart.',
     +        restart_time
         CALL WRITE_RESTART(kpp_3d_fields,kpp_const_fields,
     +        restart_outfile)
         WRITE(nuout,*) 'KPP : Wrote restart file ',restart_outfile
         CALL KPP_TIMER_TIME(kpp_timer,'Writing restart file',0)
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
      ENDIF

c     Stop opening and closing files at every output.  Just do it once at the beginning/end
c     of the simulation.
c     NPK 25/2/08 - R1
      IF (L_OUTPUT_INST) THEN
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
         CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
         CALL output_close(ncid_out)
         CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
      ENDIF
      IF (L_OUTPUT_MEAN) THEN
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
         CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
         CALL output_close(mean_ncid_out)
         CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
      ENDIF
      IF (L_OUTPUT_RANGE) THEN
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
         CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
         CALL output_close(min_ncid_out)
         CALL output_close(max_ncid_out)
         CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
         CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
      ENDIF
      WRITE(nuout,*) 'KPP : Closed output files'

#ifdef COUPLE
#ifdef OASIS2
      CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
      CALL KPP_TIMER_TIME(kpp_timer,'OASIS2 output',1)
      CALL coupled_out (kpp_3d_fields%X,ntime,.TRUE.)
      CALL KPP_TIMER_TIME(kpp_timer,'OASIS2 output',0)
c     For MPI, should do STOP 0 to return an exit code of 0
      CALL KPP_TIMER_PRINT(kpp_timer)
      STOP 0
#else
#ifdef OASIS3
c     Comment out the last coupling call, because it is hanging the coupled model
c     at the end of the integration, as the atmosphere does not post a matching
c     receive.  Depending on how we want to handle restarting the coupled integration,
c     we might want to reinstate this call, but in a modified form, so that the
c     final coupling fields from this integration can be read back in at the start
c     of the next integration.  NPK 20/10/09 - R3
c     CALL mpi1_oasis3_output(X(:,1,1),U(:,1,1),U(:,1,2))
c     CALL fluxes
      CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
      CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 output',1)
      WRITE(nuout,*) 'KPP : Calling mpi1_oasis3_terminate()'
      CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 output',0)
      CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 termination',1)
      CALL mpi1_oasis3_terminate()
      CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 termination',0)
      CALL KPP_TIMER_PRINT(kpp_timer)
c     For MPI, should do STOP 0 to return an exit code of 0
      STOP 0
#endif /*OASIS3*/
#endif /*OASIS2*/
#else
      CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
      CALL KPP_TIMER_PRINT(kpp_timer)
#endif /*COUPLE*/

      END

      SUBROUTINE initialize(kpp_3d_fields,kpp_const_fields,bottom_temp,
     +     VEC_mean,SCLR_mean,VEC_range,SCLR_range)

************************************************************************
*     Subroutine to initialize the model, some of the output is passed
*     through the common blocks
************************************************************************
c      USE kpp_type_mod
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
#include <landsea.com>
#include <constants.com>
#include <times.com>
#include <timocn.com>
#include <location.com>
#include <vert_pgrid.com>
#include <proc_swit.com>
#include <proc_pars.com>
#include <initialcon.com>
#include <ocn_advec.com>
#include <ocn_state.com>
#include <ocn_paras.com>
#include <ice_paras.com>
#include <flx_paras.com>
#include <flx_in.com>
#include <output.com>
#include <couple.com>
#include <sstclim.com>
#include <fcorr_in.com>
#include <sfcorr_in.com>
#include <relax_3d.com>
#include <bottomclim.com>
#include <currclim.com>

*     Input/output
      TYPE(kpp_3d_type),intent(inout) :: kpp_3d_fields
      TYPE(kpp_const_type),intent(inout) :: kpp_const_fields

*     Outputs
c      REAL,intent(out),allocatable :: VEC_mean(:,:,:),SCLR_mean(:,:),
c     +     bottom_temp(:)
      REAL,intent(out) :: VEC_mean(NPTS,NZP1,NVEC_MEAN),
     +     SCLR_mean(NPTS,NSCLR_MEAN),bottom_temp(NPTS),
     +     VEC_range(NPTS,NZP1,NVEC_RANGE,2),
     +     SCLR_range(NPTS,NSCLR_RANGE,2)

* Local Variablies, including some read in from name lists
      REAL dscale ! (neg)lambda parameter for defining the stretch
      REAL alat,alon,delta_lat,delta_lon
      INTEGER i,j,k,l,ipt,ix,iy
      CHARACTER*40 forcing_file

      NAMELIST/NAME_CONSTANTS/grav,vonk,sbc,twopi,onepi,TK0,spd,dpy,
     &     epsw,albocn,EL,SL,FL,FLSN
      NAMELIST/NAME_PROCSWIT/LKPP,LRI,LDD,LICE,
     &     LBIO,LNBFLX,LTGRID,LRHS,L_SSref
      NAMELIST/NAME_DOMAIN/DMAX,alon,alat,delta_lat,delta_lon,
     &     L_STRETCHGRID,dscale,L_REGGRID,L_VGRID_FILE,vgrid_file,
     &     L_SLAB,L_COLUMBIA_LAND,slab_depth
      NAMELIST/NAME_START/ L_INITDATA,initdata_file,L_INTERPINIT,
     &     L_RESTART,restart_infile
      NAMELIST/NAME_TIMES/ dtsec,startt,finalt,ndtocn
      NAMELIST/NAME_ADVEC/ L_ADVECT,advect_file,L_RELAX_SST,
     &     relax_sst_in,relax_sal_in,L_RELAX_CALCONLY,L_RELAX_SAL,
     +     L_RELAX_OCNT,relax_ocnt_in
      NAMELIST/NAME_PARAS/ paras_file,L_JERLOV
      NAMELIST/NAME_OUTPUT/ ndt_varout_inst,ndt_singout_inst,
     +     output_file,L_RESTARTW,restart_outfile,ndt_varout_mean,
     &     ndt_singout_mean,L_OUTPUT_MEAN,L_OUTPUT_INST,L_OUTPUT_RANGE,
     &     ndt_per_file,ndt_per_restart,ndt_varout_range,
     &     ndt_singout_range,zprof_varout_inst,zprof_varout_mean,
     &     zprof_varout_range,zprofs
      NAMELIST/NAME_FORCING/ L_FLUXDATA,forcing_file,L_FCORR_WITHZ,
     &     fcorrin_file,ndtupdfcorr,L_VARY_BOTTOM_TEMP,ndtupdbottom,
     &     bottomin_file,L_FCORR,L_UPD_FCORR,L_UPD_BOTTOM_TEMP,L_REST,
     &     L_PERIODIC_FCORR,L_PERIODIC_BOTTOM_TEMP,fcorr_period,
     &     L_SFCORR_WITHZ,sfcorrin_file,ndtupdsfcorr,L_SFCORR,
     &     L_UPD_SFCORR,L_PERIODIC_SFCORR,sfcorr_period,
     &     bottom_temp_period,sal_file,L_UPD_SAL,L_PERIODIC_SAL,
     &     sal_period,ndtupdsal,ocnt_file,L_UPD_OCNT,L_PERIODIC_OCNT,
     &     ocnt_period,ndtupdocnt,L_NO_FREEZE,L_NO_ISOTHERM,
     &     isotherm_bottom,isotherm_threshold,L_DAMP_CURR,dtuvdamp,
     &     L_INTERP_OCNT,ndt_interp_ocnt,L_INTERP_SAL,ndt_interp_sal
      NAMELIST/NAME_COUPLE/ L_COUPLE,ifirst,ilast,jfirst,jlast,
     &     L_CLIMSST,sstin_file,L_UPD_CLIMSST,ndtupdsst,L_CPLWGHT,
     &     cplwght_file,icein_file,L_CLIMICE,L_UPD_CLIMICE,ndtupdice,
     &     L_CLIM_ICE_DEPTH,L_CLIM_SNOW_ON_ICE,L_OUTKELVIN,
     &     L_COUPLE_CURRENTS,currin_file,L_CLIMCURR,L_UPD_CLIMCURR,
     &     ndtupdcurr,L_PERIODIC_CLIMICE,L_PERIODIC_CLIMSST,
     &     climsst_period,climice_period
      NAMELIST/NAME_LANDSEA/ L_LANDSEA,landsea_file
c
c     This is a bug fix for the IBM XLF compiler, which otherwise complains
c     about "incorrect characters" in the namelist.  If you are using the
c     IBM compiler, you need to pass -WF,-DXLF_OLDNAME when you compile
c     the KPP model.
c     NPK 5/6/09 - R2
c
#ifdef XLF_OLDNAME
      CALL SETRTEOPTS("namelist=old")
#endif
c
c     Open the namelist
c
      OPEN(75,FILE='3D_ocn.nml')
c
c     Initialse and read the constants name list
      spd=86400.                ! secs/day
      dpy=360.                  ! days/year
      twopi=8*atan(1.)          ! 2pi
      onepi=twopi/2.            ! pi
      grav=9.816                ! gravity
      vonk=0.4                  ! Von Karman's constant
      TK0=273.15                ! Kelvin of 0degC
      sbc=5.67e-8               ! Stefan Boltzmann Constant
      epsw=1.0                  ! cor.fac for departure of H2O from B.body
      albocn=0.06               ! albedo for seawater
      sice=4.0                  ! salinity of ice(?)
      EL=2.50e6                 ! Latent heat of evap. at 0C (or constant)
      SL=2512200.               ! Latent heat of evap for ice
      FL=334000.                ! Latent heat of fusion for ice
      FLSN=FL                   ! Latent heat of fusion for snow
      READ(75,NAME_CONSTANTS)
      WRITE(nuout,*) 'KPP : Read Namelist CONSTANTS'
c
c     Initialize and read the processes namelist
      LKPP=.TRUE.
      LRI=.TRUE.
      LDD=.FALSE.
      LICE=.FALSE.
      LBIO=.FALSE.
      LTGRID=.FALSE.
      LNBFLX=.FALSE.
      LRHS=.FALSE.
      L_SSref=.TRUE.
      READ(75,NAME_PROCSWIT)
      WRITE(nuout,*) 'KPP : Read Namelist PROCSWIT'
c
c     Initilalize and read the location name list
      DMAX=0.0
      alat=0.0
      alon=0.0
      delta_lat=2.5
      delta_lon=3.75
      dscale=0.0
      L_STRETCHGRID=.FALSE.
      L_REGGRID=.TRUE.
      L_VGRID_FILE=.FALSE.
      L_SLAB=.FALSE.
      L_COLUMBIA_LAND=.FALSE.
      slab_depth=0.0
      READ(75,NAME_DOMAIN)
      IF (DMAX .LE. 0.0) THEN
         WRITE(nuerr,*) 'KPP : You must specify a depth for the domain'
         CALL MIXED_ABORT
      ENDIF
      IF ((L_STRETCHGRID) .AND. (dscale .EQ. 0.0)) THEN
         WRITE(nuerr,*) 'KPP : You cannot have dscale=0 for stretched ',
     +        'grids'
         CALL MIXED_ABORT
      ENDIF
      IF ((L_SLAB) .and. (slab_depth .le. 0.0)) THEN
         WRITE(nuerr,*) 'KPP : You must specify a positive number for ',
     +        'the depth of the slab ocean.'
         CALL MIXED_ABORT
      ENDIF
      write(nuout,*) 'KPP : Read Namelist DOMAIN'
c
c     Initialize and read the landsea name list
      L_LANDSEA=.FALSE.
      READ(75,NAME_LANDSEA)
      WRITE(nuout,*) 'KPP : Read Namelist LANDSEA'
      IF (L_LANDSEA) THEN
         kpp_3d_fields%dlat(1)=alat
         kpp_3d_fields%dlon(1)=alon
         CALL init_landsea(kpp_3d_fields)
      ELSEIF (L_REGGRID) THEN
         DO iy=1,ny
            DO ix=1,nx
               ipt=(iy-1)*nx+ix
               kpp_3d_fields%dlat(ipt)=alat+(iy-1)*delta_lat
               kpp_3d_fields%dlon(ipt)=alon+(ix-1)*delta_lon
               kpp_3d_fields%ocdepth(ipt)=-10000.
               kpp_3d_fields%L_OCEAN(ipt)=.TRUE.
            ENDDO
         ENDDO
      ELSEIF (.NOT. L_REGGRID .AND. .NOT. L_LANDSEA) THEN
         WRITE(nuerr,*) 'KPP : If you set L_REGGRID=.FALSE., you must',
     +        ' specify a land-sea mask file from which to read',
     +        ' the locations of the gridpoints in the horizontal.'
      ENDIF
c
c     If coupling to the GFS, also read in the global land/sea mask on the GFS grid.
c     This allows KPP to get the global latitudes and longitudes, which it needs
c     to create a global GRIB file of SSTs to give back to the GFS.
c     NPK June 2009 - R2
c
#ifdef COUPLE
#ifdef CFS
      IF (L_LANDSEA) CALL read_landsea_global
#endif
#endif
c
c     Initialize the vertical grid
      CALL init_env(L_STRETCHGRID,dscale,kpp_3d_fields,kpp_const_fields)
      WRITE(6,*) 'after init_env zm = ',kpp_const_fields%zm
c
c     Initialize and read the start name list
      L_INITDATA= .TRUE.
      L_INTERPINIT= .TRUE.
      L_RESTART= .FALSE.
      WRITE(restart_infile,*) 'fort.30'
      READ(75,NAME_START)
      write(nuout,*) 'KPP : Read Namelist START'
c
c     Initialize and read the times namelist
      ndtocn=1
      dtsec=0.0
      startt=-999.999
      finalt=-999.999
      READ(75,NAME_TIMES)
      IF ((dtsec .LE. 0.0) .OR. (startt .LT. 0.0)
     +     .OR. (finalt .LT. 0.0)) THEN
         WRITE(nuerr,*) 'KPP : You must specify values of ',
     +        'dtsec,startt,finalt in the namelist'
         CALL MIXED_ABORT
      ENDIF
      kpp_const_fields%spd=spd
      kpp_const_fields%dtsec=dtsec
      kpp_const_fields%startt=startt*kpp_const_fields%spd
      kpp_const_fields%finalt=finalt*kpp_const_fields%spd
      kpp_const_fields%dto=kpp_const_fields%dtsec/float(ndtocn)
      nend=int((kpp_const_fields%finalt-kpp_const_fields%startt)/
     +     kpp_const_fields%dtsec)
      nstart=nint(kpp_const_fields%startt)/kpp_const_fields%dto
      IF (float(nend*ndtocn) .NE. (kpp_const_fields%finalt-
     +     kpp_const_fields%startt)/kpp_const_fields%dto) THEN
         WRITE(nuerr,*) 'KPP : The integration length is not ',
     +        'a multiple of the ocean timestep'
         WRITE(nuerr,*) 'dto=',kpp_const_fields%dto
         WRITE(nuerr,*) 'finalt=',kpp_const_fields%finalt
         WRITE(nuerr,*) 'startt=',kpp_const_fields%startt
         CALL MIXED_ABORT
      ENDIF
      kpp_const_fields%startt=kpp_const_fields%startt/
     +     kpp_const_fields%spd
      kpp_const_fields%finalt=kpp_const_fields%finalt/
     +     kpp_const_fields%spd
      kpp_const_fields%time=kpp_const_fields%startt
      WRITE(nuout,*) 'KPP : Read Namelist TIMES'
*Initialize and read the couple namelist
#ifdef COUPLE
      L_COUPLE=.TRUE.
#else
      L_COUPLE=.FALSE.
#endif
      L_COUPLE_CURRENTS=.FALSE.
      L_OUTKELVIN=.FALSE.
      L_UPD_CLIMSST=.FALSE.
      L_UPD_CLIMICE=.FALSE.
      L_CLIMICE=.FALSE.
      L_CLIMSST=.FALSE.
      L_CLIMCURR=.FALSE.
      L_BAD_ICE_DEPTH=.FALSE.
      ifirst=1
      ilast=nx
      jfirst=1
      jfirst=ny
      READ(75,NAME_COUPLE)
      write(nuout,*) 'KPP : Read Namelist COUPLE'
c      IF (L_CLIMSST) CALL read_sstin
c      IF (L_CLIMICE) CALL read_icein
c      IF (L_CLIMCURR) CALL read_surface_currents
c
c     If the model is coupled or if coupling weights
c     have been explicitly enabled, initialize the weights.
c     NPK 10/9/07 - R1
c     NPK 2/11/09 - Added #ifdef - R3
c
#ifdef COUPLE
      CALL init_cplwght(kpp_3d_fields)
#else
      IF (L_CPLWGHT) CALL init_cplwght(kpp_3d_fields)
#endif
c     Initialize and read the advection namelist
      L_ADVECT=.FALSE.
      L_RELAX_SST=.FALSE.
      L_RELAX_CALCONLY=.FALSE.
      DO iy=1,ny
         relax_sst_in(iy)=0.0
         relax_sal_in(iy)=0.0
      ENDDO
      READ(75,NAME_ADVEC)
      IF (L_ADVECT) THEN
         CALL init_advect(kpp_3d_fields)
      ELSE
         DO ipt=1,npts
            kpp_3d_fields%nmodeadv(ipt,1)=0
            kpp_3d_fields%nmodeadv(ipt,2)=0
         ENDDO
         write(nuout,*) 'KPP : No advection has been specified'
      ENDIF
      write(nuout,*) 'KPP : Read Namelist ADVEC'
      IF (L_RELAX_SST .OR. L_RELAX_SAL .OR. L_RELAX_OCNT) THEN
         CALL init_relax(kpp_3d_fields,kpp_const_fields)
      ENDIF
c     Initialize and read the paras namelist
      paras_file='3D_ocnparas.nc'
      L_JERLOV=.TRUE.
      READ(75,NAME_PARAS)
      CALL init_paras(kpp_3d_fields)
      write(nuout,*) 'KPP : Read Namelist PARAS'

c     Initialize and read the forcing namelist
      L_FLUXDATA=.FALSE.
      L_FCORR_WITHZ=.FALSE.
      L_FCORR=.FALSE.
      L_UPD_FCORR=.FALSE.
      L_SFCORR_WITHZ=.FALSE.
      L_SFCORR=.FALSE.
      L_UPD_SFCORR=.FALSE.
      L_UPD_SAL=.FALSE.
      L_VARY_BOTTOM_TEMP=.FALSE.
      L_UPD_BOTTOM_TEMP=.FALSE.
      L_REST=.FALSE.
      L_NO_FREEZE=.FALSE.
      L_NO_ISOTHERM=.FALSE.
      L_DAMP_CURR=.FALSE.
      forcing_file='1D_ocean_forcing.nc'
      ocnT_file='none'
      READ(75,NAME_FORCING)
      write(nuout,*) 'KPP : Read Namelist FORCING'
      WRITE(6,*) 'L_REST=',L_REST,'after namelist'
      IF (L_FCORR_WITHZ .AND. L_FCORR) THEN
         WRITE(nuerr,*) 'KPP : L_FCORR and L_FCORR_WITHZ are '
     &        //'mutually exclusive.  Choose one or neither.'
         CALL MIXED_ABORT
      ENDIF
      IF (L_SFCORR_WITHZ .AND. L_SFCORR) THEN
         WRITE(nuerr,*) 'KPP : L_SFCORR and L_SFCORR_WITHZ are '
     &        //'mutually exclusive.  Choose one or neither.'
         CALL MIXED_ABORT
      ENDIF
      IF (L_FCORR_WITHZ .AND. L_RELAX_SST) THEN
         WRITE(nuerr,*) 'KPP : L_FCORR_WITHZ and L_RELAX_SST are '
     &        //'mutually exclusive.  Choose one or neither.'
         CALL MIXED_ABORT
      ENDIF
c      IF (L_SFCORR_WITHZ .AND. L_RELAX_SAL) THEN
c         WRITE(nuerr,*) 'KPP : L_SFCORR_WITHZ and L_RELAX_SAL are '
c     &        //'mutually exclusive.  Choose one or neither.'
c         CALL MIXED_ABORT
c      ENDIF
      IF (L_NO_ISOTHERM .AND. (ocnT_file .eq. 'none' .or.
     &     sal_file .eq. 'none')) THEN
         WRITE(nuerr,*) 'KPP : If you specify L_NO_ISOTHERM for '
     &        //'reseting of isothermal points, you must specify files '
     &        //'from which to read climatological ocean temperature '
     &        //'(ocnT_file) and salinity (sal_file).'
         CALL MIXED_ABORT
      ELSEIF (L_NO_ISOTHERM) THEN
         kpp_const_fields%iso_bot=isotherm_bottom
         kpp_const_fields%iso_thresh=isotherm_threshold
      ENDIF
      IF (L_DAMP_CURR) THEN
         kpp_const_fields%dt_uvdamp=dtuvdamp
      ELSE
         kpp_const_fields%dt_uvdamp=10000
      ENDIF
      WRITE(6,*) kpp_3d_fields%dlon(1)
      IF (L_CLIMSST) CALL read_sstin(kpp_3d_fields,kpp_const_fields)
      IF (L_CLIMICE) CALL read_icein(kpp_3d_fields,kpp_const_fields)
      IF (L_CLIMCURR)
     +     CALL read_surface_currents(kpp_3d_fields,kpp_const_fields)
      IF (L_FCORR_WITHZ)
     +     CALL read_fcorrwithz(kpp_3d_fields,kpp_const_fields)
      IF (L_FCORR) CALL read_fcorr(kpp_3d_fields,kpp_const_fields)
      IF (L_SFCORR_WITHZ)
     +     CALL read_sfcorrwithz(kpp_3d_fields,kpp_const_fields)
      IF (L_SFCORR) CALL read_sfcorr(kpp_3d_fields,kpp_const_fields)
      IF (L_VARY_BOTTOM_TEMP)
     +     CALL read_bottom_temp(kpp_3d_fields,kpp_const_fields,
     +     bottom_temp)
      IF (L_RESTART) THEN
         CALL READ_RESTART(kpp_3d_fields,kpp_const_fields,
     +        restart_infile)
      ELSE
         CALL init_flds(kpp_3d_fields,kpp_const_fields)
         write(nuout,*) 'KPP : Temperature, salinity and currents ',
     +        ' have been initialized.'
         IF (L_UPD_BOTTOM_TEMP)
     +        CALL upd_bottom_temp(kpp_3d_fields,kpp_const_fields,
     +        bottom_temp)
      ENDIF
      CALL init_flx(kpp_3d_fields)
      IF (L_RELAX_SAL)
     +     CALL read_salinity(kpp_3d_fields,kpp_const_fields)
      IF (L_RELAX_OCNT)
     +     CALL read_ocean_temperatures(kpp_3d_fields,kpp_const_fields)
      IF (L_NO_ISOTHERM .AND. .NOT. L_RELAX_SAL
     +     .AND. .NOT. L_RELAX_OCNT) THEN
         CALL read_ocean_temperatures(kpp_3d_fields,kpp_const_fields)
         CALL read_salinity(kpp_3d_fields,kpp_const_fields)
      ENDIF
c     Currently, L_INTERP_OCNT implies L_PERIODIC_OCNT to deal with times
c     before the first time in the input file.
      IF (L_INTERP_OCNT) L_PERIODIC_OCNT=.TRUE.
      IF (L_INTERP_SAL) L_PERIODIC_SAL=.TRUE.
c
c     We need to initialize the forcing file (for atmospheric fluxes)
c     only if KPP is not coupled to an atmospheric model.
c     NPK June 2009 - R2
c
#ifndef COUPLE
      IF (L_FLUXDATA) THEN
         CALL init_flxdata(forcing_file)
      ENDIF
#endif
c     Initialize and read the output name list
      ndt_varout_inst(:)=0
      ndt_varout_mean(:)=0
      ndt_varout_range(:)=0
      zprof_varout_inst(:)=0
      zprof_varout_mean(:)=0
      zprof_varout_range(:)=0
      ndt_singout_inst(:)=0
      ndt_singout_mean(:)=0
      ndt_singout_range(:)=0
      zprofs(:,:)=0

      L_OUTPUT_MEAN=.FALSE.
      L_OUTPUT_INST=.TRUE.
      L_RESTARTW=.TRUE.
      ndt_per_restart=nend*ndtocn
      ntout_vec_inst(:)=1
      ntout_sing_inst(:)=1
      ntout_vec_mean(:)=1
      ntout_sing_mean(:)=1
      ntout_vec_range(:)=1
      ntout_sing_range(:)=1
      output_file='KPPocean'
      mean_output_file='KPPocean'
      min_output_file='KPPocean'
      max_output_file='KPPocean'
c
c     Set up defaults for ndt_per_file (timesteps between creating
c     new output files) depending on whether and how KPP is coupled.
c     GFS defaults to one day because the "coupler" operates via
c     GRIB files that need to be written/read at each coupling timestep.
c     NPK June 2009 - R2
c
#ifdef COUPLE
#ifdef CFS
      ndt_per_file=INT(kpp_const_fields%spd/(kpp_const_fields%dtsec/
     +     ndtocn))
#else
      ndt_per_file=5*INT(kpp_const_fields%spd/(kpp_const_fields%dtsec/
     +     ndtocn))
#endif /*CFS*/
#else
      ndt_per_file=5*INT(kpp_const_fields%spd/(kpp_const_fields%dtsec/
     +     ndtocn))
#endif /*COUPLE*/
      READ(75,NAME_OUTPUT)
      write(nuout,*) 'Read Namelist OUTPUT'
      zprofs_mask(:,0)=.TRUE.
      zprofs_Mask(:,1:N_ZPROFS_MAX)=.FALSE.
      zprofs_nvalid(0)=NZP1
      DO i=1,N_ZPROFS_MAX
         j=1
         DO WHILE (zprofs(j,i) .ne. 0 .and. j .le. NZP1)
            zprofs_mask(zprofs(j,i),i)=.TRUE.
            j=j+1
         END DO
         zprofs_nvalid(i)=j-1
      END DO
c
c     Set up the first output files for means and instantaneous
c     fields.
c
      flen=INDEX(output_file,' ')-1
      day_out=int(kpp_const_fields%startt+(kpp_const_fields%dtsec/
     +     ndtocn)*ndt_per_file/spd)
      write(output_file(flen+1:flen+1),'(a)') '_'
      write(output_file(flen+2:flen+6),'(i5.5)') day_out
      write(output_file(flen+7:flen+9),'(3A)') '.nc'
c
      dtout=kpp_const_fields%dto/kpp_const_fields%spd
      IF (L_OUTPUT_INST) THEN
         CALL init_output(output_file,ncid_out,kpp_3d_fields,
     +        kpp_const_fields,ndt_varout_inst,ndt_singout_inst,
     +        varid_vec,varid_sing,zprof_varout_inst,
     +        .FALSE.,.TRUE.)
         CALL output_open(output_file,ncid_out)
      ENDIF
c
      IF (L_OUTPUT_MEAN) THEN
         flen=INDEX(mean_output_file,' ')-1
         write(mean_output_file(flen+1:flen+1),'(a)') '_'
         write(mean_output_file(flen+2:flen+6),'(i5.5)') day_out
         write(mean_output_file(flen+7:flen+15),'(9A)') '_means.nc'
         WRITE(nuout,*) 'KPP : Calling init_output for '
     +        //mean_output_file
         CALL init_output(mean_output_file,mean_ncid_out,
     +        kpp_3d_fields,kpp_const_fields,ndt_varout_mean,
     +        ndt_singout_mean,varid_vec_mean,varid_sing_mean,
     +        zprof_varout_mean,.TRUE.,.FALSE.)
         CALL output_open(mean_output_file,mean_ncid_out)
      ENDIF

      IF (L_OUTPUT_RANGE) THEN
         flen=INDEX(min_output_file,' ')-1
         write(min_output_file(flen+1:flen+1),'(a)') '_'
         write(min_output_file(flen+2:flen+6),'(i5.5)') day_out
         write(min_output_file(flen+7:flen+13),'(7A)') '_min.nc'
         WRITE(nuout,*) 'KPP : Calling init_output for '
     +        //min_output_file
         CALL init_output(min_output_file,min_ncid_out,
     +        kpp_3d_fields,kpp_const_fields,ndt_varout_range,
     +        ndt_singout_range,varid_vec_range,varid_sing_range,
     +        zprof_varout_range,.FALSE.,.FALSE.)
         write(max_output_file(flen+1:flen+1),'(a)') '_'
         write(max_output_file(flen+2:flen+6),'(i5.5)') day_out
         write(max_output_file(flen+7:flen+13),'(7A)') '_max.nc'
           WRITE(nuout,*) 'KPP : Calling init_output for '
     +        //max_output_file
         CALL init_output(max_output_file,max_ncid_out,
     +        kpp_3d_fields,kpp_const_fields,ndt_varout_range,
     +        ndt_singout_range,varid_vec_range,varid_sing_range,
     +        zprof_varout_range,.FALSE.,.FALSE.)
         CALL output_open(min_output_file,min_ncid_out)
         CALL output_open(max_output_file,max_ncid_out)
      ENDIF
c
c     Call routine to copy constants and logicals needed for ocean
c     physics into the kpp_const_fields derived type.  Added for
c     compatability with OpenMP DEFAULT(private). NPK 8/2/13
      CALL kpp_const_fields_init(kpp_const_fields)
c
      kpp_const_fields%ntime=0
      CALL wm_ws_lookup(kpp_const_fields)
      CALL init_ocn(kpp_3d_fields,kpp_const_fields)
c     Write out the data from the initial condition
      IF ( .NOT. L_RESTART .AND. L_OUTPUT_INST) THEN
         DO l=1,N_VAROUTS
            IF (ndt_varout_inst(l) .gt. 0)
     +           CALL output_inst(kpp_3d_fields,kpp_const_fields,
     +           l,varid_vec(l),zprof_varout_inst(l),ntout_vec_inst(l))
         ENDDO
         DO l=1,N_SINGOUTS
            IF (ndt_singout_inst(l) .gt. 0)
     +           CALL output_inst(kpp_3d_fields,kpp_const_fields,
     +           l+N_VAROUTS,varid_sing(l),0,ntout_sing_inst(l))
         ENDDO
      ENDIF

c     Set the means to zero initially
      VEC_mean(:,:,:) = 0.
      SCLR_mean(:,:) = 0.
c     Set ranges to large values
      VEC_range(:,:,:,1)=2E20
      SCLR_range(:,:,1)=2E20
      VEC_range(:,:,:,2)=-2E20
      SCLR_range(:,:,2)=-2E20

      CLOSE(75)
      WRITE(6,*) 'Returning from initialise'

      RETURN
      END

      SUBROUTINE WRITE_RESTART(kpp_3d_fields,kpp_const_fields,
     +     restart_outfile)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
c#include <times.com>
c#include <ocn_paras.com>
c#include <ocn_state.com>
c#include <output.com>
c#include <flx_in.com>
c
c     Inputs
c
      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      CHARACTER(LEN=17) :: restart_outfile

c     When the number of points in the model (NX*NY*NZP1) becomes
c     quite large, we exceed the maximum size for Fortran unformatted
c     binary files on certain machines.  The IF test below
c     works around this by splitting the restart file in two.
c     %Us and %Xs are the largest fields, so they get their own file.

      WRITE(6,*) 'Total number of points = ',REAL(NPTS)*REAL(NZP1)
      IF ( REAL(NPTS)*REAL(NZP1) .LT. 3000000. ) THEN
         OPEN(31,FILE=restart_outfile,status='unknown',
     +        form='unformatted')
         WRITE(31) kpp_const_fields%time,kpp_3d_fields%U,
     +        kpp_3d_fields%X,kpp_3d_fields%CP,
     +        kpp_3d_fields%rho,kpp_3d_fields%hmix,kpp_3d_fields%kmix,
     +        kpp_3d_fields%Sref,kpp_3d_fields%SSref,
     +        kpp_3d_fields%Ssurf,
     +        kpp_3d_fields%Tref,kpp_3d_fields%old,kpp_3d_fields%new,
     +        kpp_3d_fields%Us,kpp_3d_fields%Xs,kpp_3d_fields%hmixd
         CLOSE(31)
      ELSE
         OPEN(31,FILE=restart_outfile//'.1',status='unknown',
     +        form='unformatted')
         WRITE(31) kpp_const_fields%time,kpp_3d_fields%U,
     +        kpp_3d_fields%X,kpp_3d_fields%CP,
     +        kpp_3d_fields%rho,kpp_3d_fields%hmix,kpp_3d_fields%kmix,
     +        kpp_3d_fields%Sref,kpp_3d_fields%SSref,
     +        kpp_3d_fields%Ssurf,
     +        kpp_3d_fields%Tref,kpp_3d_fields%old,kpp_3d_fields%new
         CLOSE(31)
         OPEN(32,FILE=restart_outfile//'.2',status='unknown',
     +        form='unformatted')
         WRITE(32) kpp_3d_fields%Us,kpp_3d_fields%Xs,kpp_3d_fields%hmixd
         CLOSE(32)
      ENDIF

      RETURN
      END

      SUBROUTINE READ_RESTART(kpp_3d_fields,kpp_const_fields,
     +     restart_infile)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
c
c     Inputs
c
      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      CHARACTER(LEN=17) :: restart_infile

      WRITE(6,*) 'Total number of points = ',REAL(NPTS)*REAL(NZP1)
      IF ( REAL(NPTS)*REAL(NZP1) .LT. 3000000. ) THEN
         OPEN(30,FILE=restart_infile,status='unknown',
     +        form='unformatted')
         READ(30) kpp_const_fields%time,kpp_3d_fields%U,
     +        kpp_3d_fields%X,kpp_3d_fields%CP,
     +        kpp_3d_fields%rho,kpp_3d_fields%hmix,kpp_3d_fields%kmix,
     +        kpp_3d_fields%Sref,kpp_3d_fields%SSref,
     +        kpp_3d_fields%Ssurf,
     +        kpp_3d_fields%Tref,kpp_3d_fields%old,kpp_3d_fields%new,
     +        kpp_3d_fields%Us,kpp_3d_fields%Xs,kpp_3d_fields%hmixd
         CLOSE(30)
      ELSE
         OPEN(30,FILE=restart_infile//'.1',status='unknown',
     +        form='unformatted')
         READ(30) kpp_const_fields%time,kpp_3d_fields%U,
     +        kpp_3d_fields%X,kpp_3d_fields%CP,
     +        kpp_3d_fields%rho,kpp_3d_fields%hmix,kpp_3d_fields%kmix,
     +        kpp_3d_fields%Sref,kpp_3d_fields%SSref,
     +        kpp_3d_fields%Ssurf,
     +        kpp_3d_fields%Tref,kpp_3d_fields%old,kpp_3d_fields%new
         CLOSE(30)
         OPEN(31,FILE=restart_infile//'.2',status='unknown',
     +        form='unformatted')
         READ(31) kpp_3d_fields%Us,kpp_3d_fields%Xs,kpp_3d_fields%hmixd
         CLOSE(31)
      ENDIF

      IF (abs(kpp_const_fields%time-kpp_const_fields%startt) .GT. 1.e-4)
     +     THEN
         WRITE(nuerr,*) 'Start time doesn''t match the restart record'
         WRITE(nuerr,*) 'Start time in restart record = ',
     +        kpp_const_fields%time
         WRITE(nuerr,*) 'Start time in namelist = ',
     +        kpp_const_fields%startt
c         CALL MIXED_ABORT
	 kpp_const_fields%time=kpp_const_fields%startt
      ENDIF

      RETURN
      END

      SUBROUTINE MIXED_ABORT

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

c Close output files so we are not left with
c un-readable netCDF output. NPK 18/5/13
c      CALL output_close(ncid_out)
c      CALL mean_output_close(mean_ncid_out)

#ifdef COUPLE
#ifdef OASIS2
c
c     Support for stopping in OASIS2
c     NPK April 2009 - R2
c
      CALL halte('KPP : MIXED_ABORT called')
#else
c
c     Support for stopping in OASIS3
c     NPK 2/11/09 - R3
c
#ifdef OASIS3
      CALL mpi1_oasis3_terminate()
#else
#ifdef CFS
c
c     Support for stopping with the CFS coupler
c     Unsure how to stop the model for the GFS - Just stop?
c     NPK June 2009 - R2
c
      STOP
#endif /*CFS*/
#endif /*OASIS3*/
#endif /*OASIS2*/
#else
      STOP
#endif /*COUPLE*/

      END

      SUBROUTINE init_relax(kpp_3d_fields,kpp_const_fields)
c
c     Re-write logic to allow for relaxing either SST or
c     salinity - NPK 24/08/11
c
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
      include 'constants.com'
      include 'ocn_advec.com'
      include 'couple.com'
      include 'sstclim.com'
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
               kpp_3d_fields%relax_sst(ipoint)=1./(relax_sst_in(iy)*
     +              kpp_const_fields%spd)
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
               kpp_3d_fields%relax_sal(ipoint)=1./(relax_sal_in(iy)*
     +              kpp_const_fields%spd)
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
               kpp_3d_fields%relax_ocnT(ipoint)=1./(relax_ocnT_in(iy)*
     +              kpp_const_fields%spd)
            ENDDO
         ENDIF
      ENDDO
      CALL upd_sst0(kpp_3d_fields)
      DO iy=1,ny
         DO ix=1,nx
            ipoint=(iy-1)*nx+ix
            kpp_3d_fields%fcorr(ipoint)=0.0
            kpp_3d_fields%scorr(ipoint,:)=0.0
         ENDDO
      ENDDO

      write(nuout,*) 'calculated SST0, fcorr and scorr'

      RETURN
      END

      SUBROUTINE upd_sst0(kpp_3d_fields)

c     Written by NPK 27/8/07

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
      include 'constants.com'
      include 'ocn_advec.com'
      include 'couple.com'
      include 'sstclim.com'

      INTEGER ix,iy,ipoint

      TYPE(kpp_3d_type) :: kpp_3d_fields
      REAL sst_in(NX_GLOBE,NY_GLOBE,1)
      COMMON /save_sstin/ sst_in

      DO iy=1,ny
         DO ix=1,nx
            ipoint=(iy-1)*nx+ix
c            SST0(ipoint)=SST_in(ix+ifirst_sst-1,iy+jfirst_sst-1,1)
            kpp_3d_fields%SST0(ipoint)=SST_in(ix+ifirst-1,iy+jfirst-1,1)
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE upd_bottom_temp(kpp_3d_fields,kpp_const_fields,
     +     bottom_temp)

c     Written by NPK 10/4/08

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <kpp_3d_type.com>
#include <landsea.com>

      INTEGER ipt,z
      REAL bottom_temp(NPTS)
      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

      DO ipt=1,npts
         kpp_3d_fields%tinc_fcorr(ipt,NZP1)=
     +        bottom_temp(ipt)-kpp_3D_fields%X(ipt,NZP1,1)
         kpp_3d_fields%ocnTcorr(ipt,NZP1)=
     +        kpp_3d_fields%tinc_fcorr(ipt,NZP1)*
     +        kpp_3d_fields%rho(ipt,NZP1)*kpp_3d_fields%cp(ipt,NZP1)/
     +        kpp_const_fields%dto
         kpp_3d_fields%X(ipt,NZP1,1) = bottom_temp(ipt)
      ENDDO

      RETURN
      END

      SUBROUTINE check_freezing(kpp_3d_fields)
c
c     Check whether the temperature at any (x,z) point is less than the
c     threshold for sea ice (-1.8C).  If it is, reset it to -1.8C and
c     set a flag.  The flag can be requested as a diagnostic (singout 9).
c     Note that the value of the flag is equal to the *fraction* of levels
c     at that point that were < -1.8C.
c
      IMPLICIT NONE
#include <kpp_3d_type.com>
      INTEGER ipt,z

      TYPE(kpp_3d_type) :: kpp_3d_fields

      DO ipt=1,npts
         IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
            DO z=1,NZP1
               IF (kpp_3d_fields%X(ipt,z,1) .lt. -1.8) THEN
                  kpp_3d_fields%tinc_fcorr(ipt,z)=
     +                 kpp_3d_fields%tinc_fcorr(ipt,z)+
     +                 (-1.8-kpp_3d_fields%X(ipt,z,1))
                  kpp_3d_fields%X(ipt,z,1)=-1.8
                  kpp_3d_fields%freeze_flag(ipt)=
     +                 kpp_3d_fields%freeze_flag(ipt)+1.0/REAL(NZP1)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      RETURN
      END

      SUBROUTINE check_isothermal(kpp_3d_fields,kpp_const_fields)
c
c     Check whether the temperature difference between the surface
c     and a user-specified level (presumably deep) is less than a user-specified
c     threshold (presumably small).  If so, reset the temperature and salinity
c     profiles to climatological values.  Added to prevent spurious very
c     deep mixing that creates unrealistic isothermal (and isohaline) layers.
c
c     NPK 15/5/2013 for R4.
c
      IMPLICIT NONE
#include <kpp_3d_type.com>
      INTEGER ipt,z,j
      REAL dz_total,dtdz_total,dz

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

#ifdef OPENMP
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(kpp_3d_fields,kpp_const_fields)
!$OMP DO SCHEDULE(dynamic)
#endif
      DO ipt=1,npts
         IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
            dtdz_total=0.
            dz_total=0.
            DO j=2,kpp_const_fields%iso_bot
               dz=kpp_const_fields%zm(j)-kpp_const_fields%zm(j-1)
               dtdz_total=dtdz_total+
     +              ABS((kpp_3d_fields%X(ipt,j,1)-
     +              kpp_3d_fields%X(ipt,j-1,1)))*dz
               dz_total=dz_total+dz
            ENDDO
            dtdz_total=dtdz_total/dz_total
c If resetting to climatology because of isothermal layer (rather than because of
c computational instability trap in ocn.f), then set reset_flag to a negative
c value (-1*number of interations in of semi-implicit integration in ocn.f).
            IF (ABS(dtdz_total).lt.kpp_const_fields%iso_thresh) THEN
               kpp_3d_fields%X(ipt,:,1)=kpp_3d_fields%ocnT_clim(ipt,:)
               kpp_3d_fields%X(ipt,:,2)=kpp_3d_fields%sal_clim(ipt,:)
               kpp_3d_fields%reset_flag(ipt)=(-1.)*
     +              kpp_3d_fields%reset_flag(ipt)
            ENDIF
         ELSE
            kpp_3d_fields%reset_flag(ipt)=0
         ENDIF
      ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif

      RETURN
      END

      SUBROUTINE interp_ocnT(kpp_3d_fields,kpp_const_fields)
      IMPLICIT NONE
#include <kpp_3d_type.com>
#include <relax_3d.com>
      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      INTEGER prev_time,next_time,true_time
      REAL prev_weight,next_weight,ndays_upd_ocnT
      REAL, allocatable :: prev_ocnT(:,:),next_ocnT(:,:)

      allocate(prev_ocnT(NPTS,NZP1))
      allocate(next_ocnT(NPTS,NZP1))
      true_time=kpp_const_fields%time
      ndays_upd_ocnT=ndtupdocnT*kpp_const_fields%dto/
     +     kpp_const_fields%spd

!     Read ocean temperatures for previous time
      prev_time=FLOOR((true_time+ndays_upd_ocnT/2)/ndays_upd_ocnT)*
     +     ndays_upd_ocnT-ndays_upd_ocnT*0.5
      IF (prev_time .lt. 0) THEN
         prev_weight=(ndays_upd_ocnT-ABS(true_time-prev_time))/
     +        ndays_upd_ocnT
         prev_time=prev_time+ocnT_period
      ELSE
         prev_weight=(ndays_upd_ocnT-(true_time-prev_time))/
     +        ndays_upd_ocnT
      ENDIF
      WRITE(6,*) 'interp_ocnT : true_time = ',true_time
      WRITE(6,*) 'interp_ocnT : prev_time = ',prev_time
      WRITE(6,*) 'interp_ocnT : prev_weight = ',prev_weight
      kpp_const_fields%time=prev_time
      CALL READ_OCEAN_TEMPERATURES(kpp_3d_fields,kpp_const_fields)
      prev_ocnT=kpp_3d_fields%ocnT_clim

!     Read ocean temperatures for next time
      next_time=prev_time+ndays_upd_ocnT
      next_weight=1-prev_weight
      WRITE(6,*) 'interp_ocnT : next_time = ',next_time
      WRITE(6,*) 'interp_ocnT : next_weight = ',next_weight
      kpp_const_fields%time=next_time
      CALL READ_OCEAN_TEMPERATURES(kpp_3d_fields,kpp_const_fields)
      next_ocnT=kpp_3d_fields%ocnT_clim

      kpp_3d_fields%ocnT_clim=next_ocnT*next_weight+
     +     prev_ocnT*prev_weight
      kpp_const_fields%time=true_time

      RETURN
      END

      SUBROUTINE interp_sal(kpp_3d_fields,kpp_const_fields)
      IMPLICIT NONE
#include <kpp_3d_type.com>
#include <relax_3d.com>
      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      INTEGER prev_time,next_time,true_time
      REAL prev_weight,next_weight,ndays_upd_sal
      REAL, allocatable :: prev_sal(:,:),next_sal(:,:)

      allocate(prev_sal(NPTS,NZP1))
      allocate(next_sal(NPTS,NZP1))
      true_time=kpp_const_fields%time
      ndays_upd_sal=ndtupdsal*kpp_const_fields%dto/
     +     kpp_const_fields%spd

!     Read ocean salinity for previous time
      prev_time=FLOOR((true_time+ndays_upd_sal/2)/ndays_upd_sal)*
     +     ndays_upd_sal-ndays_upd_sal*0.5
      IF (prev_time .lt. 0) THEN
         prev_weight=(ndays_upd_sal-ABS(true_time-prev_time))/
     +        ndays_upd_sal
         prev_time=prev_time+sal_period
      ELSE
         prev_weight=(ndays_upd_sal-(true_time-prev_time))/
     +        ndays_upd_sal
      ENDIF
      WRITE(6,*) 'interp_sal : true_time = ',true_time
      WRITE(6,*) 'interp_sal : prev_time = ',prev_time
      WRITE(6,*) 'interp_sal : prev_weight = ',prev_weight
      kpp_const_fields%time=prev_time
      CALL READ_SALINITY(kpp_3d_fields,kpp_const_fields)
      prev_sal=kpp_3d_fields%sal_clim

!     Read ocean salinity for next time
      next_time=prev_time+ndays_upd_sal
      next_weight=1-prev_weight
      WRITE(6,*) 'interp_sal : next_time = ',next_time
      WRITE(6,*) 'interp_sal : next_weight = ',next_weight
      kpp_const_fields%time=next_time
      CALL READ_SALINITY(kpp_3d_fields,kpp_const_fields)
      next_sal=kpp_3d_fields%sal_clim

      kpp_3d_fields%sal_clim=next_sal*next_weight+
     +     prev_sal*prev_weight
      kpp_const_fields%time=true_time

      RETURN
      END
