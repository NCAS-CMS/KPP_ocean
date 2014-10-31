PROGRAM mckpp_driver

#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type),allocatable :: kpp_3d_fields
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  TYPE(kpp_timer_type) :: kpp_timer
  
  INTEGER :: ntime
  
  allocate(kpp_3d_fields)
  
  ! Only master processor does this
  ! Initialize model simulation (constants, initial and boundary conditions, namelist options)
  WRITE(6,*) 'MCKPP_DRIVER: Calling MCKPP_INITIALIZE'
  CALL MCKPP_INITIALIZE(kpp_3d_fields,kpp_const_fields)
  WRITE(6,*) 'MCKPP_DRIVER: Returned from MCKPP_INITIALIZE'
  ! Need to broadcast wm_ws lookup table
  
#ifdef COUPLE
  ! Initialize coupling to atmosphere
  CALL MCKPP_COUPLING_INITIALIZE(kpp_3d_fields,kpp_const_fields)
#endif
  
  ! Main time-stepping loop
  DO ntime=1,kpp_const_fields%nend*kpp_const_fields%ndtocn
     kpp_const_fields%ntime=ntime
     
     ! Update surface fluxes if necessary, by invoking coupler if coupled, or from file if forced
     IF (MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndtocn) .EQ. 0) THEN
        WRITE(6,*) 'MCKPP_DRIVER: Calling MCKPP_FLUXES_UPDATE'
        CALL MCKPP_FLUXES_UPDATE(kpp_3d_fields,kpp_const_fields,kpp_timer)
        WRITE(6,*) 'MCKPP_DRIVER: Returned from MCKPP_FLUXES_UPDATE'
     ENDIF
          
     ! Update other (non-coupler-related) boundary conditions if necessary, from files
     CALL MCKPP_BOUNDARY_UPDATE(kpp_3d_fields,kpp_const_fields,kpp_timer)
     
     ! Ocean physics and associated checks
     CALL MCKPP_PHYSICS_DRIVER(kpp_3d_fields,kpp_const_fields,kpp_timer)
     
     ! Checkpointing
     CALL MCKPP_RESTART_CONTROL(kpp_3d_fields,kpp_const_fields,kpp_timer)
     
     ! Output management
     CALL MCKPP_OUTPUT_CONTROL(kpp_3d_fields,kpp_const_fields,kpp_timer)
     
     ! If coupled, return ocean fields to atmosphere on coupling timestep
#ifdef COUPLE
     IF ((MOD(ntime,kpp_const_fields%ndtocn) .EQ. 0) .AND. ntime .NE. nend*ndtocn) THEN
        CALL MCKPP_COUPLING_STEP(kpp_3d_fields,kpp_const_fields)
     ENDIF
#endif /*COUPLE*/             
  ENDDO
  
  ! Terminate model at end of integration
  CALL MCKPP_FINALIZE(kpp_3d_fields,kpp_const_fields,kpp_timer)
  
END PROGRAM mckpp_driver
