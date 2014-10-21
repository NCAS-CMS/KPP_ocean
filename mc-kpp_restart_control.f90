SUBROUTINE mckpp_restart_control(kpp_3d_fields,kpp_const_fields,kpp_timer)
  IMPLICIT NONE

#include <mc-kpp_3d_type.com>
  
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  TYPE(kpp_timer_type) :: kpp_timer

  CHARACTER*5 :: restart_time
  CHARACTER*17 :: restart_outfile
  
  IF (kpp_const_fields%L_RESTARTW) THEN
     ! Write restart file either every ndt_per_restart timesteps
     ! or at the end of the simulation.
     IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndt_per_restart) .EQ. 0 .OR. &
          kpp_const_fields%ndt_per_restart .EQ. kpp_const_fields%nend*kpp_const_fields%ndtocn) THEN
        !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
        !CALL KPP_TIMER_TIME(kpp_timer,'Writing Restart File',1)
        IF (kpp_const_fields%time .lt. 10) THEN
           WRITE(restart_time,'(A4,I1)') '0000',FLOOR(kpp_const_fields%time)
        ELSEIF (kpp_const_fields%time .lt. 100) THEN
           WRITE(restart_time,'(A3,I2)') '000',FLOOR(kpp_const_fields%time)
        ELSEIF (kpp_const_fields%time .lt. 1000) THEN
           WRITE(restart_time,'(A2,I3)') '00',FLOOR(kpp_const_fields%time)
        ELSEIF (kpp_const_fields%time .lt. 10000) THEN
           WRITE(restart_time,'(A1,I4)') '0',FLOOR(kpp_const_fields%time)
        ELSE
           WRITE(restart_time,'(I5)') FLOOR(kpp_const_fields%time)
        ENDIF
        WRITE(restart_outfile,'(A12,A5)') 'KPP.restart.',restart_time
        CALL MCKPP_RESTART_IO_WRITE(kpp_3d_fields,kpp_const_fields,restart_outfile)
        !CALL KPP_TIMER_TIME(kpp_timer,'Writing Restart File',0)
        !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
     ENDIF
  ENDIF
  
  RETURN
END SUBROUTINE mckpp_restart_control
