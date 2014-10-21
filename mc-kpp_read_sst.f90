SUBROUTINE mckpp_read_sst(kpp_3d_fields,kpp_const_fields)
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)

#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

#ifdef COUPLE
  integer,parameter :: sst_nx=NX_GLOBE,sst_ny=NY_GLOBE
#else
  integer,parameter :: sst_nx=NX,sst_ny=NY
#endif
  REAL :: offset_sst
  INTEGER :: status,ncid
  REAL*4 :: var_in(sst_nx,sst_ny,1),time_in,first_timein,last_timein,sstclim_time,&
       longitudes(NX_GLOBE),latitudes(NY_GLOBE)
  INTEGER :: varid, time_varid,lat_varid,lon_varid,lon_dimid,lat_dimid,time_dimid
  INTEGER count(3),start(3)
  INTEGER ix,iy,nlat_file,nlon_file,ntime_file
  CHARACTER(LEN=30) tmp_name
  
  status=NF_OPEN(kpp_const_fields%sst_file,0,ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)

  count(1)=sst_nx
  count(2)=sst_ny
  count(3)=1
  start(1)=1
  start(2)=1
  start(3)=1
  
  status=NF_INQ_VARID(ncid,'sst',varid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
#ifndef COUPLE
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(ncid,'SST climatology','latitude','longitude','t',kpp_3d_fields%dlon(1),&
       kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
#else
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(ncid,'t',time_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIMID(ncid,'t',time_dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIM(ncid,time_dimid,tmp_name,ntime_file)
  status=NF_GET_VAR1_REAL(ncid,time_varid,start(3),first_timein)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VAR1_REAL(ncid,time_varid,ntime_file,last_timein)
#endif
  sstclim_time=kpp_const_fields%time+0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdsst
  IF (sstclim_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_CLIMSST) THEN 
        DO WHILE (sstclim_time .gt. last_timein)
           sstclim_time=sstclim_time-kpp_const_fields%climsst_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'Time for which to read SST exceeds the last time in the netCDF file &
             & and L_PERIODIC_CLIMSST has not been specified.  Attempting to read SST will lead to &
             & an error, so aborting now ...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF
  write(nuout,*) 'Reading climatological SST for time ',sstclim_time
  
  start(3)=NINT((sstclim_time-first_timein)*kpp_const_fields%spd/&
       (kpp_const_fields%dto*kpp_const_fields%ndtupdsst))+1

  status=NF_GET_VAR1_REAL(ncid,time_varid,start(3),time_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  IF (abs(time_in-sstclim_time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     write(nuerr,*) 'Cannot find time,',sstclim_time,'in SST climatology file'
     write(nuerr,*) 'The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
 
  status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)  
  status=NF_CLOSE(ncid)

  ! KPP expects temperatures in CELSIUS.  If climatological SSTs are
  ! in Kelvin, subtract 273.15.
  offset_sst=0.
  DO ix=1,sst_nx
     DO iy=1,sst_ny
        IF (var_in(ix,iy,1) .gt. 200 .and. var_in(ix,iy,1) .lt. 400) &
             offset_sst = 273.15
     END DO
  END DO
  
  DO ix=1,sst_nx
     DO iy=1,sst_ny
        kpp_3d_fields%sst(ix,iy)=var_in(ix,iy,1)-offset_sst          
        IF (.NOT. kpp_const_fields%L_CLIMICE) kpp_3d_fields%iceconc(ix,iy)=0.0
        IF (.NOT. kpp_const_fields%L_CLIMCURR) THEN
           kpp_3d_fields%usf(ix,iy)=0.0
           kpp_3d_fields%vsf(ix,iy)=0.0
        ENDIF
     ENDDO
  ENDDO
  
  WRITE(nuout,*) 'Finished read_sstin'
  
END SUBROUTINE mckpp_read_sst
