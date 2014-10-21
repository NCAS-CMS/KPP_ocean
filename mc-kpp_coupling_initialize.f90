SUBROUTINE mckpp_coupling_initialize(kpp_3d_fields,kpp_const_fields)
  IMPLICIT NONE

#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  INTEGER,parameter :: nuout=6, nuerr=0

#ifdef OASIS2
#include <param.h>
#endif

#ifdef OASIS3
#include <kpp_oasis3.inc>
#endif

#ifdef OASIS2
  IF ((im .NE. NX_GLOBE) .OR. (jm .NE. NY_GLOBE)) THEN
     
     ! Test that KPP and OASIS have the same horizontal grid sizes.
     ! These are controlled by NX_GLOBE and NY_GLOBE in parameter.inc
     ! (for KPP) and by im and jm in param.h (for OASIS).
     ! This applies to OASIS2 only.
     WRITE(nuout,*) 'KPP : KPP and OASIS2 do not have the same horizontal grid sizes.  KPP has NX_GLOBE=',&
          NX_GLOBE,'and NY_GLOBE=',NY_GLOBE,'while OASIS2 has im=',im,'and jm=',jm,'.  You can control the KPP settings in ',&
          'parameter.inc and the OASIS2 settings in param.h.'
     CALL halte('im,jm not equal to NX_GLOBE,NY_GLOBE')
  ELSE
     ! Initialize the OASIS2 coupling interface
     !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
     !CALL KPP_TIMER_TIME(kpp_timer,'OASIS2 initialization',1)
     CALL inicmo(kpp_const_fields%nend*kpp_const_fields%ndtocn,kpp_const_fields%ndtocn,int(kpp_const_fields%dto)) 
     !CALL KPP_TIMER_TIME(kpp_timer,'OASIS2 initialization',0)
     !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
  ENDIF
#endif /*OASIS2*/

#ifdef OASIS3
  CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
  CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 initialization',1)
  CALL MCKPP_COUPLING_OASIS3_INITIALIZE()
  CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 initialization',0)
  CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)

  ! When coupling to the UM via OASIS3, we must send the initial ocean fields to
  ! the atmosphere.  The UM does an OASIS3 "get" first, then an OASIS3 "send."
  ! We must match this with a "send" before we post our first "get", otherwise
  ! the whole system deadlocks and no one has any fun.  
  
  ! Note that the OASIS3 "get" is handled by the first call to <fluxes> in
  ! the time-stepping loop below.     NPK 2/10/09
  CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
  CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 output',1)
  CALL MCKPP_COUPLING_OASIS3_OUTPUT(kpp_3d_fields,kpp_const_fields)
  CALL KPP_TIMER_TIME(kpp_timer,'OASIS3 output',0)
  CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
#endif /*OASIS3*/

  RETURN
END SUBROUTINE mckpp_coupling_initialize
