SUBROUTINE mckpp_fluxes_update(kpp_3d_fields,kpp_const_fields,kpp_timer)
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
  ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
!#include <landsea.com>
#ifdef COUPLE
#ifdef OASIS3
#include <kpp_oasis3.inc>
#endif /*OASIS3*/
#endif /*COUPLE*/
  
!  include 'ocn_paras.com'
!  include 'flx_paras.com'
!  include 'flx_in.com'
!  include 'local_pt.com'
!#include <timocn.com>
!#include <initialcon.com>
!      include 'times.com'
!      include 'couple.com'

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  TYPE(kpp_timer_type) :: kpp_timer
  
  !      REAL taux(NPTS),tauy(NPTS),
  !     $     swf(NPTS),lwf(NPTS),lhf(NPTS),shf(NPTS),
  !     $     rain(NPTS),snow(NPTS)
  CHARACTER(LEN=19) :: trans_timer_name
  INTEGER :: ipt
#ifdef OPENMP
  INTEGER tid,OMP_GET_THREAD_NUM      
#endif

#ifdef COUPLE
#ifdef OASIS2
  call coupled_flux(kpp_3d_fields%swf,kpp_3d_fields%lwf,kpp_3d_fields%rain,kpp_const_fields%ntime-1)
  call coupled_stress(kpp_3d_fields%taux,kpp_3d_fields%tauy,kpp_const_fields%ntime-1)
#else
#ifdef CFS
  CALL read_gfs_forcing(kpp_3d_fields%swf,kpp_3d_fields%lwf,kpp_3d_fields%rain,kpp_3d_fields%taux,kpp_3d_fields%tauy)
#else
#ifdef OASIS3
  ! Normal coupling - no writing to or reading from netCDF files
  CALL MCKPP_COUPLING_OASIS3_INPUT(kpp_3d_fields,kpp_const_fields)
  
  ! HadGEM3 passes zeros at the first timestep for a new run (i.e., NRUN)
  ! Thus, if this is NOT a restart run, we need to provide a file
  ! of fluxes for KPP for the first coupling timestep.      
  ! I am specifying "kpp_initfluxes.nc" as the filename here, although
  ! we could always add a namelist option later to control this.
  ! NPK 15/10/09, revised 6/11/09 to specify .NOT. L_RESTART
  ! as HadGEM3 does pass good fields for a restart run (i.e., CRUN)
  IF (kpp_const_fields%ntime .EQ. 1 .AND. .NOT. kpp_const_fields%L_RESTART) THEN 
     CALL MCKPP_INITIALIZE_FLUXES_FILE('kpp_initfluxes.nc')
     CALL MCKPP_READ_FLUXES(kpp_3d_fields,kpp_const_fields)
     ! Convert to variables expected for a coupled model
     DO ipt=1,npts
        IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
           kpp_3d_fields%lwf(ipt)=kpp_3d_fields%lwf(ipt)+kpp_3d_fields%lhf(ipt)+&
                kpp_3d_fields%shf(ipt)-kpp_3d_fields%snow(ipt)*kpp_const_fields%FLSN
        ENDIF
     ENDDO
  ENDIF
#endif /*OASIS3*/
#endif /*CFS*/
#endif /*OASIS2*/
! All coupled models do this step      
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(shared) 
!$OMP& PRIVATE(ipt,kpp_1d_fields,tid,trans_timer_name)
  tid=OMP_GET_THREAD_NUM()
  WRITE(trans_timer_name,'(A17,I2)') 'KPP 3D/2D thread ',tid
#else
  WRITE(trans_timer_name,'(A19)') 'KPP 3D/2D thread 01'
#endif
  !CALL kpp_timer_time(kpp_timer,trans_timer_name,1)
#ifdef OPENMP
!$OMP DO SCHEDULE(dynamic)
#endif
  DO ipt=1,npts
     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
        IF ((kpp_3d_fields%taux(ipt) .EQ. 0.0) .AND. (kpp_3d_fields%tauy(ipt) .EQ. 0.0)) &
           kpp_3d_fields%taux(ipt)=1.e-10
        kpp_3d_fields%sflux(ipt,1,5,0)=kpp_3d_fields%taux(ipt)
        kpp_3d_fields%sflux(ipt,2,5,0)=kpp_3d_fields%tauy(ipt)
        kpp_3d_fields%sflux(ipt,3,5,0)=kpp_3d_fields%swf(ipt)               
        kpp_3d_fields%sflux(ipt,4,5,0)=kpp_3d_fields%lwf(ipt) 
        kpp_3d_fields%sflux(ipt,5,5,0)=0.0 ! Melting of sea-ice = 0.0               
        kpp_3d_fields%sflux(ipt,6,5,0)=kpp_3d_fields%rain(ipt) ! assuming rain = P-E
        CALL mckpp_fields_3dto1d(kpp_3d_fields,ipt,kpp_1d_fields) 
        call mckpp_fluxes_ntflx(kpp_1d_fields,kpp_const_fields)
        CALL mckpp_fields_1dto3d(kpp_1d_fields,ipt,kpp_3d_fields)
     ENDIF
  ENDDO
#ifdef OPENMP
!$OMP END DO
#endif
  !CALL kpp_timer_time(kpp_timer,trans_timer_name,0)
#ifdef OPENMP
!$OMP END PARALLEL
#endif
#else  /* NOT COUPLED */
  ! Get fluxes for the forced case
  IF (.NOT. kpp_const_fields%L_FLUXDATA) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(shared) PRIVATE(ipt)
!$OMP DO SCHEDULE(static)
#endif
     DO ipt=1,npts
        kpp_3d_fields%taux(ipt)=0.01
        kpp_3d_fields%tauy(ipt)=0.0
        kpp_3d_fields%swf(ipt)=200.0
        kpp_3d_fields%lwf(ipt)=0.0
        kpp_3d_fields%lhf(ipt)=-150.0
        kpp_3d_fields%shf(ipt)=0.0
        kpp_3d_fields%rain(ipt)=6e-5
        kpp_3d_fields%snow(ipt)=0.0
     ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  ELSE
     call MCKPP_READ_FLUXES(kpp_3d_fields,kpp_const_fields)
  ENDIF
      
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(shared) PRIVATE(ipt,kpp_1d_fields)
!$OMP DO SCHEDULE(dynamic)
#endif
  DO ipt=1,npts         
     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN 
        IF ((kpp_3d_fields%taux(ipt) .EQ. 0.0) .AND. (kpp_3d_fields%tauy(ipt) .EQ. 0.0)) &
             kpp_3d_fields%taux(ipt)=1.e-10
        IF (.NOT. kpp_const_fields%L_REST) THEN
           kpp_3d_fields%sflux(ipt,1,5,0)=kpp_3d_fields%taux(ipt)
           kpp_3d_fields%sflux(ipt,2,5,0)=kpp_3d_fields%tauy(ipt)
           kpp_3d_fields%sflux(ipt,3,5,0)=kpp_3d_fields%swf(ipt)            
           kpp_3d_fields%sflux(ipt,4,5,0)=kpp_3d_fields%lwf(ipt)+&
                kpp_3d_fields%lhf(ipt)+kpp_3d_fields%shf(ipt)-&
                kpp_3d_fields%snow(ipt)*kpp_const_fields%FLSN
           kpp_3d_fields%sflux(ipt,6,5,0)=(kpp_3d_fields%rain(ipt)+&
                kpp_3d_fields%snow(ipt)+(kpp_3d_fields%lhf(ipt)/kpp_const_fields%EL))
           kpp_3d_fields%sflux(ipt,5,5,0)=1e-10 ! Melting of sea-ice = 0.0
        ELSE
           kpp_3d_fields%sflux(ipt,1,5,0)=1.e-10
           kpp_3d_fields%sflux(ipt,2,5,0)=0.00
           kpp_3d_fields%sflux(ipt,3,5,0)=300.00
           kpp_3d_fields%sflux(ipt,4,5,0)=-300.00
           kpp_3d_fields%sflux(ipt,5,5,0)=0.00
           kpp_3d_fields%sflux(ipt,6,5,0)=0.00
        ENDIF
        CALL mckpp_fields_3dto1d(kpp_3d_fields,ipt,kpp_1d_fields)
        call mckpp_fluxes_ntflux(kpp_1d_fields,kpp_const_fields)
        CALL mckpp_fields_1dto3d(kpp_1d_fields,ipt,kpp_3d_fields)
     ENDIF
  ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
#endif /*COUPLE*/
  RETURN
END SUBROUTINE mckpp_fluxes_update
