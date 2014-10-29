SUBROUTINE MCKPP_READ_FCORR_2D(kpp_3d_fields,kpp_const_fields)
  IMPLICIT NONE
  INTEGER nuout,nuerr,start(3),count(3)
  INTEGER ix,iy,ipoint,fcorr_varid,status,lat_varid,lon_varid,time_varid,&
       lat_dimid,lon_dimid,time_dimid,fcorr_ncid,k,nlat_file,nlon_file,ntime_file
  PARAMETER (nuout=6,nuerr=0)
     
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  REAL*4 ixx,jyy,fcorr_twod_in(NX,NY,1),latitudes(NY),longitudes(NX),z(NZP1),&
       first_timein,time_in,fcorr_time,ndays_upd_fcorr,last_timein
  CHARACTER(LEN=30) tmp_name
  
  ! Read in a NetCDF file containing a time-varying flux correction
  ! at the surface only.  Frequency of read is controlled by ndtupdfcorr
  ! in the namelist
  ! NPK 29/06/08

  status = NF_OPEN(kpp_const_fields%fcorr_file,0,fcorr_ncid)
  IF (status.NE.0) CALL MCKPP_HANDLE_ERR(status)
  WRITE(nuout,*) 'Opened flux-correction input file'
  
  count(1)=NX
  count(2)=NY      
  count(3)=1
  start(1)=1
  start(2)=1
  start(3)=1
  
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(fcorr_ncid,'flux correction','latitude','longitude',&
       't',kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
  
  status=NF_INQ_VARID(fcorr_ncid,'fcorr',fcorr_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
       
  ndays_upd_fcorr = kpp_const_fields%ndtupdfcorr*kpp_const_fields%dto/kpp_const_fields%spd
  fcorr_time=(ndays_upd_fcorr)*(FLOOR(kpp_const_fields%time)*NINT(kpp_const_fields%spd)/&
       (kpp_const_fields%ndtupdfcorr*NINT(kpp_const_fields%dto)))+&
       (0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdfcorr)
  
  IF (fcorr_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_FCORR) THEN 
        DO WHILE (fcorr_time .gt. last_timein)
           fcorr_time=fcorr_time-kpp_const_fields%fcorr_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'Time for which to read the flux corrections exceeds the last time in the netCDF file &
             & and L_PERIODIC_FCORR has not been specified. Attempting to read flux corrections will lead to &
             & an error, so aborting now ...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF

  write(nuout,*) 'Reading flux correction for time ',fcorr_time
  start(3)=NINT((fcorr_time-first_timein)*kpp_const_fields%spd/(kpp_const_fields%dto*kpp_const_fields%ndtupdfcorr))+1
  write(nuout,*) 'Flux corrections are being read from position',start(3)
  status=NF_GET_VAR1_REAL(fcorr_ncid,time_varid,start(3),time_in)
  
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-fcorr_time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     write(nuerr,*) 'Cannot find time',fcorr_time,'in flux-correction input file'
     write(nuerr,*) 'The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
  status=NF_GET_VARA_REAL(fcorr_ncid,fcorr_varid,start,count,fcorr_twod_in)
  write(nuout,*) 'Flux corrections have been read from position',start(3)
  
  !    Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
  !    into one long array with dimension NPTS.         
  DO ix=1,NX
     DO iy=1,NY
        ipoint=(iy-1)*nx+ix
        kpp_3d_fields%fcorr_twod(ipoint)=fcorr_twod_in(ix,iy,1)
     ENDDO
  ENDDO
  
  status=NF_CLOSE(fcorr_ncid)
  
END SUBROUTINE MCKPP_READ_FCORR_2D

SUBROUTINE MCKPP_READ_FCORR_3D(kpp_3d_fields,kpp_const_fields)
  
  IMPLICIT NONE
  INTEGER nuout,nuerr,start(4),count(4)
  INTEGER ix,iy,iz,ipoint,fcorr_varid,status,lat_varid,lon_varid,z_varid,z_dimid,time_varid,&
       fcorr_ncid,k,lat_dimid,lon_dimid,time_dimid,nlon_file,nlat_file,ntime_file,nz_file
  
  PARAMETER (nuout=6,nuerr=0)
#include <netcdf.inc>
  ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  REAL*4 ixx,jyy,first_timein,time_in,fcorr_time,ndays_upd_fcorr,last_timein
  CHARACTER(LEN=30) tmp_name
  REAL*4, allocatable :: fcorr_in(:,:,:,:),longitudes(:),latitudes(:),z(:)
  
  ! Read in a NetCDF file containing a 
  ! time-varying flux correction at every model vertical level.
  ! Frequency of read is controlled by ndtupdfcorr in the namelist
  ! NPK 12/02/08
  allocate(fcorr_in(NX,NY,NZP1,1))
  allocate(longitudes(NX_GLOBE))
  allocate(latitudes(NY_GLOBE))
  allocate(z(NZP1))
  
  WRITE(nuout,*) 'MCKPP_READ_FCORR_3D: Opening flux-correction file ',kpp_const_fields%fcorr_file
  status=NF_OPEN(kpp_const_fields%fcorr_file,0,fcorr_ncid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  WRITE(nuout,*) 'Opened flux-correction input file ',fcorr_ncid
  
  count=(/NX,NY,NZP1,1/)
  start=(/1,1,1,1/)
  
  status=NF_INQ_VARID(fcorr_ncid,'z',z_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)  
  status=NF_INQ_DIMID(fcorr_ncid,'z',z_dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIM(fcorr_ncid,z_dimid,tmp_name,nz_file)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (NZP1.ne.nz_file) THEN
     WRITE(nuout,*) 'Input file for flux corrections does not have the correct number of vertical levels. ',&
          'It should have ',NZP1,' but instead has ',nz_file
     CALL MCKPP_ABORT
  ELSE
     status=NF_GET_VAR_REAL(fcorr_ncid,z_varid,z)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     WRITE(nuout,*) 'Read in depths from the flux-correction input file'
  ENDIF
  
  WRITE(6,*) 'MCKPP_READ_FCORR_3D: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(fcorr_ncid,'flux correction','latitude','longitude',&
       't',kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
  status=NF_INQ_VARID(fcorr_ncid,'fcorr',fcorr_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  ndays_upd_fcorr = kpp_const_fields%ndtupdfcorr*kpp_const_fields%dto/kpp_const_fields%spd
  WRITE(nuout,*) ndays_upd_fcorr,FLOOR(kpp_const_fields%time)*NINT(kpp_const_fields%spd),&
       kpp_const_fields%ndtupdfcorr*NINT(kpp_const_fields%dto),&
       0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdfcorr
  fcorr_time=(ndays_upd_fcorr)*FLOOR(kpp_const_fields%time)*NINT(kpp_const_fields%spd)/&
       FLOAT(kpp_const_fields%ndtupdfcorr*NINT(kpp_const_fields%dto))+&
       (0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdfcorr)
  WRITE(nuout,*) fcorr_time,last_timein
      
  IF (fcorr_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_FCORR) THEN 
        DO WHILE (fcorr_time .gt. last_timein)
           fcorr_time=fcorr_time-kpp_const_fields%fcorr_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'Time for which to read the flux corrections exceeds the &
             & last time in the netCDF file and L_PERIODIC_FCORR has not been specified. &
             & Attempting to read flux corrections will lead to an error, so aborting now ...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF
  
  write(nuout,*) 'Reading flux correction for time ',fcorr_time
  start(4)=NINT((fcorr_time-first_timein)*kpp_const_fields%spd/&
       (kpp_const_fields%dto*kpp_const_fields%ndtupdfcorr))+1
  write(nuout,*) 'Flux corrections are being read from position',start(4)
  status=NF_GET_VAR1_REAL(fcorr_ncid,time_varid,start(4),time_in)
      
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-fcorr_time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     write(nuerr,*) 'Cannot find time',fcorr_time,'in flux-correction input file'
     write(nuerr,*) 'The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
  status=NF_GET_VARA_REAL(fcorr_ncid,fcorr_varid,start,count,fcorr_in)
  write(nuout,*) 'Flux corrections have been read from position',start(4)      
  
!     Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
!     into one long array with dimension NPTS.         
  DO ix=1,NX
     DO iy=1,NY
        ipoint=(iy-1)*nx+ix
        DO k=1,NZP1
           kpp_3d_fields%fcorr_withz(ipoint,k)=fcorr_in(ix,iy,k,1)
        ENDDO
     ENDDO
  ENDDO
  
  status=NF_CLOSE(fcorr_ncid)
  deallocate(fcorr_in)
  deallocate(longitudes)
  deallocate(latitudes)
  deallocate(z)
  
  RETURN
END SUBROUTINE MCKPP_READ_FCORR_3D
