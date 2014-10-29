SUBROUTINE mckpp_physics_driver(kpp_3d_fields,kpp_const_fields,kpp_timer)
  IMPLICIT NONE

#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  TYPE(kpp_timer_type) :: kpp_timer
  
  TYPE(kpp_1d_type) :: kpp_1d_fields
  INTEGER, parameter :: nuout=6, nuerr=0
  INTEGER :: ipt
#ifdef OPENMP
  INTEGER :: tid,OMP_GET_THREAD_NUM
#endif
  CHARACTER(LEN=21) phys_timer_name
  CHARACTER(LEN=19) trans_timer_name
  
  kpp_const_fields%time=kpp_const_fields%startt+kpp_const_fields%ntime*&
       kpp_const_fields%dto/kpp_const_fields%spd
  
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(kpp_3d_fields,kpp_const_fields) &
!$OMP SHARED(kpp_timer) & 
!$OMP PRIVATE(trans_timer_name,phys_timer_name,tid)
  tid=OMP_GET_THREAD_NUM()
  WRITE(trans_timer_name,'(A17,I2)') 'KPP 3D/1D thread ',tid
  WRITE(phys_timer_name,'(A19,I2)') 'KPP Physics thread ',tid
!$OMP DO SCHEDULE(dynamic)
#else
  WRITE(trans_timer_name,'(A19)') 'KPP 3D/1D thread 01'
  WRITE(phys_timer_name,'(A21)') 'KPP Physics thread 01'
#endif
  DO ipt=1,npts
     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
        !CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,1)
        CALL mckpp_fields_3dto1d(kpp_3d_fields,ipt,kpp_1d_fields)
        !CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,0)
        !CALL KPP_TIMER_TIME(kpp_timer,phys_timer_name,1)
                
        CALL mckpp_physics_ocnstep(kpp_1d_fields,kpp_const_fields)       
        CALL mckpp_physics_overrides_check_profile(kpp_1d_fields,kpp_const_fields)        
        
        !CALL KPP_TIMER_TIME(kpp_timer,phys_timer_name,0)
        !CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,1)
        CALL mckpp_fields_1dto3d(kpp_1d_fields,ipt,kpp_3d_fields)
        !CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,0)
     ENDIF
  ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  
  IF (kpp_const_fields%L_VARY_BOTTOM_TEMP) THEN
     !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
     !CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
     CALL MCKPP_PHYSICS_OVERRIDES_BOTTOMTEMP(kpp_3d_fields,kpp_const_fields)
     !CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
     !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
  ENDIF

  RETURN
END SUBROUTINE mckpp_physics_driver
