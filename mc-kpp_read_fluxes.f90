SUBROUTINE mckpp_read_fluxes(kpp_3d_fields,kpp_const_fields)
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  REAL*4 time_in,time
      
  INTEGER dimid,varid
  INTEGER ipt,ix,iy,nx_in,ny_in
  REAL*4, allocatable :: x_in(:),y_in(:),var_in(:,:)
      
  INTEGER status,start(3),count(3)
  
  allocate(x_in(NX_GLOBE))
  allocate(y_in(NY_GLOBE))
  allocate(var_in(NX,NY))
  
  start(1)=1
  start(2)=1
  start(3)=1
  count(1)=1
  count(2)=1
  count(3)=1
  
  status=NF_INQ_DIMID(kpp_const_fields%flx_ncid,'longitude',dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIMLEN (kpp_const_fields%flx_ncid,dimid,nx_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'longitude',varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VAR_REAL(kpp_const_fields%flx_ncid,varid,x_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ix=1
  DO WHILE (abs(x_in(ix)-kpp_3d_fields%dlon(1)) .GT. 1.e-3)
     ix=ix+1
     IF (ix .GE. nx_in) THEN
        write(nuerr,*) 'Error reading fluxes'
        write(nuerr,*) 'Can''t find longitude ',kpp_3d_fields%dlon(1),&
             ' in range ',x_in(1),x_in(nx_in)
        CALL MCKPP_ABORT
     ENDIF
  ENDDO
  start(1)=ix
  count(1)=nx
  
  status=NF_INQ_DIMID(kpp_const_fields%flx_ncid,'latitude',dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIMLEN (kpp_const_fields%flx_ncid,dimid,ny_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'latitude',varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VAR_REAL(kpp_const_fields%flx_ncid,varid,y_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  iy=1
  DO WHILE (abs(y_in(iy)-kpp_3d_fields%dlat(1)) .GT. 1.e-3)
     iy=iy+1
     IF (iy .GE. ny_in) THEN
        write(nuerr,*) 'Error reading fluxes'
        write(nuerr,*) 'Can''t find latitude ',kpp_3d_fields%dlat(1),&
             ' in range ',y_in(1),y_in(ny_in)
        CALL MCKPP_ABORT
     ENDIF
  ENDDO
  start(2)=iy
  count(2)=ny
      
  time=kpp_const_fields%time+0.5*kpp_const_fields%dtsec/kpp_const_fields%spd
  write(nuout,*) ' Reading fluxes for time ',time,kpp_const_fields%time,&
       kpp_const_fields%dtsec,kpp_const_fields%spd
  start(3)=MAX(NINT((time-kpp_const_fields%flx_first_timein)*kpp_const_fields%spd/kpp_const_fields%dtsec)+1,1)
  WRITE(6,*) 'Reading time from time point ',start(3)
  WRITE(6,*) 'first_timein=',kpp_const_fields%flx_first_timein,'time=',time
  status=NF_GET_VAR1_REAL(kpp_const_fields%flx_ncid,kpp_const_fields%flx_timein_id,start(3),time_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#ifndef OASIS3      
  IF (abs(time_in-time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     write(nuerr,*) 'Cannot find time,',time,'in fluxes file'
     write(nuerr,*) 'The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
#endif
  WRITE(6,*) 'Reading fluxes from time point ',start(3)
  status=NF_GET_VARA_REAL(kpp_const_fields%flx_ncid,kpp_const_fields%flx_varin_id(1),start,count,var_in)
  WRITE(6,*) 'Read taux'
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
     DO ix=1,nx
        ipt=(iy-1)*nx+ix            
        kpp_3d_fields%taux(ipt)=var_in(ix,iy)
     ENDDO
  ENDDO
  status=NF_GET_VARA_REAL(kpp_const_fields%flx_ncid,kpp_const_fields%flx_varin_id(2),start,count,var_in)
  WRITE(6,*) 'Read tauy'
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
     DO ix=1,nx
        ipt=(iy-1)*nx+ix
        kpp_3d_fields%tauy(ipt)=var_in(ix,iy)
     ENDDO
  ENDDO
  status=NF_GET_VARA_REAL(kpp_const_fields%flx_ncid,kpp_const_fields%flx_varin_id(3),start,count,var_in)
  WRITE(6,*) 'Read swf'
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
     DO ix=1,nx
        ipt=(iy-1)*nx+ix
        kpp_3d_fields%swf(ipt)=var_in(ix,iy)
     ENDDO
  ENDDO
  status=NF_GET_VARA_REAL(kpp_const_fields%flx_ncid,kpp_const_fields%flx_varin_id(4),start,count,var_in)
  WRITE(6,*) 'Read lwf'
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
     DO ix=1,nx
        ipt=(iy-1)*nx+ix
        kpp_3d_fields%lwf(ipt)=var_in(ix,iy)
     ENDDO
  ENDDO
  status=NF_GET_VARA_REAL(kpp_const_fields%flx_ncid,kpp_const_fields%flx_varin_id(5),start,count,var_in)
  WRITE(6,*) 'Read lhf'
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
     DO ix=1,nx
        ipt=(iy-1)*nx+ix
        kpp_3d_fields%lhf(ipt)=var_in(ix,iy)
     ENDDO
  ENDDO
  status=NF_GET_VARA_REAL(kpp_const_fields%flx_ncid,kpp_const_fields%flx_varin_id(6),start,count,var_in)
  WRITE(6,*) 'Read shf'
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)      
  DO iy=1,ny
     DO ix=1,nx
        ipt=(iy-1)*nx+ix
        kpp_3d_fields%shf(ipt)=var_in(ix,iy)
     ENDDO
  ENDDO
  status=NF_GET_VARA_REAL(kpp_const_fields%flx_ncid,kpp_const_fields%flx_varin_id(7),start,count,var_in)
  WRITE(6,*) 'Read rain'
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
     DO ix=1,nx
        ipt=(iy-1)*nx+ix
        kpp_3d_fields%rain(ipt)=var_in(ix,iy)
     ENDDO
  ENDDO
  DO iy=1,ny
     DO ix=1,nx
        ipt=(iy-1)*nx+ix
        kpp_3d_fields%snow(ipt)=0.0
     ENDDO
  ENDDO
  WRITE(6,*) 'Set snow to zero'
  WRITE(6,*) 'Finished reading fluxes'
  
  RETURN
END SUBROUTINE mckpp_read_fluxes
