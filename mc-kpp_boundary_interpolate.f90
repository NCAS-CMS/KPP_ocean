SUBROUTINE MCKPP_BOUNDARY_INTERPOLATE_TEMP(kpp_3d_fields,kpp_const_fields)      
  IMPLICIT NONE

#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER prev_time,next_time,true_time
  REAL prev_weight,next_weight,ndays_upd_ocnT
  REAL, allocatable :: prev_ocnT(:,:),next_ocnT(:,:)
  
  allocate(prev_ocnT(NPTS,NZP1))
  allocate(next_ocnT(NPTS,NZP1))
  true_time=kpp_const_fields%time
  ndays_upd_ocnT=kpp_const_fields%ndtupdocnT*kpp_const_fields%dto/kpp_const_fields%spd
      
  ! Read ocean temperatures for previous time
  prev_time=FLOOR((true_time+ndays_upd_ocnT/2)/ndays_upd_ocnT)*ndays_upd_ocnT-ndays_upd_ocnT*0.5
  IF (prev_time .lt. 0) THEN
     prev_weight=(ndays_upd_ocnT-ABS(true_time-prev_time))/ndays_upd_ocnT
     prev_time=prev_time+kpp_const_fields%ocnT_period
  ELSE
     prev_weight=(ndays_upd_ocnT-(true_time-prev_time))/ndays_upd_ocnT
  ENDIF
  WRITE(6,*) 'interp_ocnT : true_time = ',true_time
  WRITE(6,*) 'interp_ocnT : prev_time = ',prev_time
  WRITE(6,*) 'interp_ocnT : prev_weight = ',prev_weight
  kpp_const_fields%time=prev_time
  CALL MCKPP_READ_TEMPERATURES_3D(kpp_3d_fields,kpp_const_fields)
  prev_ocnT=kpp_3d_fields%ocnT_clim

  ! Read ocean temperatures for next time
  next_time=prev_time+ndays_upd_ocnT
  next_weight=1-prev_weight
  WRITE(6,*) 'interp_ocnT : next_time = ',next_time
  WRITE(6,*) 'interp_ocnT : next_weight = ',next_weight
  kpp_const_fields%time=next_time
  CALL MCKPP_READ_TEMPERATURES_3D(kpp_3d_fields,kpp_const_fields)
  next_ocnT=kpp_3d_fields%ocnT_clim
  
  kpp_3d_fields%ocnT_clim=next_ocnT*next_weight+prev_ocnT*prev_weight
  kpp_const_fields%time=true_time

  RETURN
END SUBROUTINE MCKPP_BOUNDARY_INTERPOLATE_TEMP

SUBROUTINE MCKPP_BOUNDARY_INTERPOLATE_SAL(kpp_3d_fields,kpp_const_fields)
  IMPLICIT NONE

#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER prev_time,next_time,true_time
  REAL prev_weight,next_weight,ndays_upd_sal
  REAL, allocatable :: prev_sal(:,:),next_sal(:,:)
  
  allocate(prev_sal(NPTS,NZP1))
  allocate(next_sal(NPTS,NZP1))
  true_time=kpp_const_fields%time
  ndays_upd_sal=kpp_const_fields%ndtupdsal*kpp_const_fields%dto/kpp_const_fields%spd
  
  ! Read ocean salinity for previous time
  prev_time=FLOOR((true_time+ndays_upd_sal/2)/ndays_upd_sal)*ndays_upd_sal-ndays_upd_sal*0.5
  IF (prev_time .lt. 0) THEN
     prev_weight=(ndays_upd_sal-ABS(true_time-prev_time))/ndays_upd_sal
     prev_time=prev_time+kpp_const_fields%sal_period
  ELSE
     prev_weight=(ndays_upd_sal-(true_time-prev_time))/ndays_upd_sal
  ENDIF
  WRITE(6,*) 'interp_sal : true_time = ',true_time
  WRITE(6,*) 'interp_sal : prev_time = ',prev_time
  WRITE(6,*) 'interp_sal : prev_weight = ',prev_weight
  kpp_const_fields%time=prev_time
  CALL MCKPP_READ_SALINITY_3D(kpp_3d_fields,kpp_const_fields)
  prev_sal=kpp_3d_fields%sal_clim
  
  ! Read ocean salinity for next time
  next_time=prev_time+ndays_upd_sal
  next_weight=1-prev_weight
  WRITE(6,*) 'interp_sal : next_time = ',next_time
  WRITE(6,*) 'interp_sal : next_weight = ',next_weight
  kpp_const_fields%time=next_time
  CALL MCKPP_READ_SALINITY_3D(kpp_3d_fields,kpp_const_fields)
  next_sal=kpp_3d_fields%sal_clim
  
  kpp_3d_fields%sal_clim=next_sal*next_weight+prev_sal*prev_weight
  kpp_const_fields%time=true_time
  
  RETURN
END SUBROUTINE MCKPP_BOUNDARY_INTERPOLATE_SAL
