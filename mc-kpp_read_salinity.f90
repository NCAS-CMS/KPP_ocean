SUBROUTINE MCKPP_READ_SALINITY_3D(kpp_3d_fields,kpp_const_fields)
  IMPLICIT NONE
  INTEGER nuout,nuerr,start(4),count(4)
  INTEGER ix,iy,iz,ipoint,sal_varid,status,lat_varid,lon_varid,z_varid,z_dimid,time_varid,&
       sal_ncid,k,lat_dimid,lon_dimid,time_dimid,nlon_file,nlat_file,ntime_file,nz_file

  PARAMETER (nuout=6,nuerr=0)
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
      
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  REAL*4 ixx,jyy,first_timein,time_in,sal_time,ndays_upd_sal,last_timein
  CHARACTER(LEN=30) tmp_name
  REAL*4, allocatable :: sal_in(:,:,:,:),latitudes(:),longitudes(:),z(:)

! Read in a NetCDF file containing a 
! time-varying salinity field at every model vertical level.
! Frequency of read is controlled by ndtupdsal in the namelist
! NPK 12/02/08
     
  allocate(sal_in(NX,NY,NZP1,1))
  allocate(longitudes(NX_GLOBE))
  allocate(latitudes(NY_GLOBE))
  allocate(z(NZP1))
  
  WRITE(nuout,*) 'Trying to open salinity input file ',kpp_const_fields%sal_file
  status=NF_OPEN(kpp_const_fields%sal_file,0,sal_ncid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  WRITE(nuout,*) 'Opened salinity input file ',sal_ncid
  
  count=(/NX,NY,NZP1,1/)
  start=(/1,1,1,1/)
  
  status=NF_INQ_VARID(sal_ncid,'z',z_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)  
  status=NF_INQ_DIMID(sal_ncid,'z',z_dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIM(sal_ncid,z_dimid,tmp_name,nz_file)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (NZP1.ne.nz_file) THEN
     WRITE(nuout,*) 'Input file for salinity climatology does not have the correct number of vertical levels. ',&
          'It should have ',NZP1,' but instead has ',nz_file
     CALL MCKPP_ABORT
  ELSE
     status=NF_GET_VAR_REAL(sal_ncid,z_varid,z)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     WRITE(nuout,*) 'Read in depths from the salinity climatology input file'
  ENDIF
  
  WRITE(6,*) kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1)
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(sal_ncid,'salinity clim','latitude','longitude','t',kpp_3d_fields%dlon(1),&
       kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
  
  status=NF_INQ_VARID(sal_ncid,'salinity',sal_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  ndays_upd_sal = kpp_const_fields%ndtupdsal*kpp_const_fields%dto/kpp_const_fields%spd
  sal_time=(ndays_upd_sal)*(FLOOR(kpp_const_fields%time)*NINT(kpp_const_fields%spd)/&
       (kpp_const_fields%ndtupdsal*NINT(kpp_const_fields%dto)))+&
       (0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdsal)
  WRITE(nuout,*) sal_time,last_timein
  
  IF (sal_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_SAL) THEN 
        DO WHILE (sal_time .gt. last_timein)
           sal_time=sal_time-kpp_const_fields%sal_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'Time for which to read the salinity climatology exceeds the last time in the netCDF file &
             & and L_PERIODIC_SAL has not been specified. Attempting to read salinity climatology will lead to &
             & an error, so aborting now ...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF
  
  write(nuout,*) 'Reading salinity for time ',sal_time
  start(4)=NINT((sal_time-first_timein)*kpp_const_fields%spd/(kpp_const_fields%dto*kpp_const_fields%ndtupdsal))+1
  write(nuout,*) 'Salinity values are being read from position',start(4)
  status=NF_GET_VAR1_REAL(sal_ncid,time_varid,start(4),time_in)
  
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-sal_time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     write(nuerr,*) 'Cannot find time',sal_time,'in salinity climatology input file'
     write(nuerr,*) 'The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
  status=NF_GET_VARA_REAL(sal_ncid,sal_varid,start,count,sal_in)
  write(nuout,*) 'Salinity climatology data have been read from position',start(4)

  ! Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
  ! into one long array with dimension NPTS.      
  DO ix=1,NX
     DO iy=1,NY
        ipoint=(iy-1)*nx+ix
        DO k=1,NZP1
           ! Subtract reference salinity from climatology, for compatability with
           ! salinity values stored in main model.  
           kpp_3d_fields%sal_clim(ipoint,k)=sal_in(ix,iy,k,1)-kpp_3d_fields%Sref(ipoint)
        ENDDO
     ENDDO
  ENDDO
  
  status=NF_CLOSE(sal_ncid)
  deallocate(sal_in)
  deallocate(longitudes)
  deallocate(latitudes)
  deallocate(z)
  
  RETURN
END SUBROUTINE MCKPP_READ_SALINITY_3D
