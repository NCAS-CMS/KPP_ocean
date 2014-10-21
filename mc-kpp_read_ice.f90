SUBROUTINE mckpp_read_ice(kpp_3d_fields,kpp_const_fields)
  
  ! Read in ice concentrations from a user-provided netCDF file.
  ! Called _only_ if L_CLIMICE is TRUE in 3D_ocn.nml.
  ! Probably only necessary when coupling the model to an atmosphere
  ! that requires realistic ice concentrations.
  ! Written by Nick Klingaman, 11/01/08.
  
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
  ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
  !  include 'constants.com'
  !  include 'couple.com'
  !  include 'times.com'
  !  include 'timocn.com'
  !  include 'sstclim.com'      
  
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
#ifdef COUPLE
  integer,parameter :: ice_nx=NX_GLOBE,ice_ny=NY_GLOBE
#else
  integer,parameter :: ice_nx=NX,ice_ny=NY
#endif     
  REAL :: max_ice,min_ice
  REAL*4 :: var_in(ice_nx,ice_ny,1),iceclim_time,first_timein,last_timein,&
       time_in,latitudes(NY_GLOBE),longitudes(NX_GLOBE)
  INTEGER count(3),start(3)
  INTEGER ix,iy,status,ncid,varid,time_varid,lon_varid,lat_varid,time_dimid,&
       lon_dimid,lat_dimid,ntime_file,nlon_file,nlat_file
  CHARACTER(LEN=30) tmp_name
      
  ! Set start and count to read a global field if coupled,
  ! or a regional field if not coupled.
  count(1)=ice_nx
  count(2)=ice_ny
  count(3)=1
  start(1)=1
  start(2)=1
  start(3)=1

  ! Open the netCDF file and find the correct latitude,
  ! longitude and time.
  WRITE(nuout,*) 'Opening ice input file'
  status=NF_OPEN(kpp_const_fields%ice_file,0,ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)

#ifndef COUPLE
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(ncid,'ice climatology','latitude','longitude',&
       't',kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
#else
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(ncid,'t',time_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIMID(ncid,'t',time_dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIM(ncid,time_dimid,tmp_name,ntime_file)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VAR1_REAL(ncid,time_varid,start(3),first_timein)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VAR1_REAL(ncid,time_varid,ntime_file,last_timein)
#endif      

  status=NF_INQ_VARID(ncid,'iceconc',varid)
  iceclim_time=kpp_const_fields%time+0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdice       
  IF (iceclim_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_CLIMICE) THEN 
        DO WHILE (iceclim_time .gt. last_timein)
           iceclim_time=iceclim_time-kpp_const_fields%climice_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'Time for which to read ice exceeds the last time in the netCDF file &
             & and L_PERIODIC_CLIMICE has not been specified.  Attempting to read ice will lead to &
             & an error, so aborting now ...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF
      
  write(nuout,*) 'Reading climatological ICECONC for time ',iceclim_time
  start(3)=NINT((iceclim_time-first_timein)*kpp_const_fields%spd/&
       (kpp_const_fields%dto*kpp_const_fields%ndtupdice))+1      
  status=NF_GET_VAR1_REAL(ncid,time_varid,start(3),time_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-iceclim_time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     write(nuerr,*) 'Cannot find time,',iceclim_time,'in ice concentration climatology file'
     write(nuerr,*) 'The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
  write(nuout,*) 'Ice concentrations are being read from position',start(3)
  WRITE(nuout,*) 'Start = ',start,'Count = ',count
  status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  WRITE(nuout,*) 'Ice concentrations have been read from position',start(3)
  
  max_ice = -1000.
  min_ice = 1000.
  DO ix=1,ice_nx
     DO iy=1,ice_ny
        kpp_3d_fields%iceconc(ix,iy) = var_in(ix,iy,1)
        IF (kpp_3d_fields%iceconc(ix,iy) .gt. max_ice) max_ice = kpp_3d_fields%iceconc(ix,iy)
        IF (kpp_3d_fields%iceconc(ix,iy) .lt. min_ice) min_ice = kpp_3d_fields%iceconc(ix,iy)
     ENDDO
  ENDDO
  
  IF (kpp_const_fields%L_CLIM_ICE_DEPTH) THEN 
     status=NF_INQ_VARID(ncid,'icedepth',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     WRITE(nuout,*) 'Reading climatological ICEDEPTH for time ',iceclim_time        
     WRITE(nuout,*) 'Ice depths are being read from position', start(3)
     WRITE(nuout,*) 'Start = ',start,'Count = ',count
     status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)      
     DO ix=1,ice_nx
        DO iy=1,ice_ny
           kpp_3d_fields%icedepth(ix,iy)=var_in(ix,iy,1)
        ENDDO
     ENDDO
  ENDIF

  IF (kpp_const_fields%L_CLIM_SNOW_ON_ICE) THEN
     status=NF_INQ_VARID(ncid,'snowdepth',varid)
     IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     WRITE(nuout,*) 'Reading climatological SNOWDEPTH for time ',iceclim_time
     WRITE(nuout,*) 'Snow depths on sea ice are being read from position', start(3)
     WRITE(nuout,*) 'Start = ',start,'Count = ',count
     status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
     IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)     
     DO ix=1,ice_nx
        DO iy=1,ice_ny
           kpp_3d_fields%snowdepth(ix,iy)=var_in(ix,iy,1)
        ENDDO
     ENDDO
  ENDIF
  
  status=NF_CLOSE(ncid)
  
  RETURN
END SUBROUTINE mckpp_read_ice
