SUBROUTINE MCKPP_READ_TEMPERATURES_3D(kpp_3d_fields,kpp_const_fields)
  IMPLICIT NONE

#include <mc-kpp_3d_type.com>
#include <netcdf.inc>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER nuout,nuerr,start(4),count(4)
  INTEGER ix,iy,iz,ipoint,ocnT_varid,status,lat_varid,lon_varid,z_varid,z_dimid,time_varid,&
       ocnT_ncid,k,lat_dimid,lon_dimid,time_dimid,nlon_file,nlat_file,ntime_file,nz_file,prev_start
  
  PARAMETER (nuout=6,nuerr=0)
  REAL*4, allocatable :: ocnT_in(:,:,:,:),latitudes(:),longitudes(:),z(:)
  REAL*4 ixx,jyy,first_timein,time_in,ocnT_time,ndays_upd_ocnT,last_timein
  CHARACTER(LEN=30) tmp_name

  allocate(ocnT_in(NX,NY,NZP1,1))
  allocate(longitudes(NX_GLOBE))
  allocate(latitudes(NY_GLOBE))
  allocate(z(NZP1))
  
  WRITE(nuout,*) 'Trying to open ocean temperature input file ',kpp_const_fields%ocnT_file
  status=NF_OPEN(kpp_const_fields%ocnT_file,0,ocnT_ncid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  WRITE(nuout,*) 'Opened ocean temperature input file ',ocnT_ncid
  
  count=(/NX,NY,NZP1,1/)
  start=(/1,1,1,1/)
  
  status=NF_INQ_VARID(ocnT_ncid,'z',z_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)  
  status=NF_INQ_DIMID(ocnT_ncid,'z',z_dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIM(ocnT_ncid,z_dimid,tmp_name,nz_file)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (NZP1.ne.nz_file) THEN
     WRITE(nuout,*) 'Input file for ocean temperature climatology does not have the &
          & correct number of vertical levels. It should have ',NZP1,' but instead has ',nz_file
     CALL MCKPP_ABORT
  ELSE
     status=NF_GET_VAR_REAL(ocnT_ncid,z_varid,z)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     WRITE(nuout,*) 'Read in depths from the ocean temperature climatology input file'
  ENDIF

  WRITE(6,*) kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1)
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(ocnT_ncid,'ocean temp clim','latitude','longitude','t',kpp_3d_fields%dlon(1),&
       kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
  
  status=NF_INQ_VARID(ocnT_ncid,'temperature',ocnT_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  ndays_upd_ocnT = kpp_const_fields%ndtupdocnT*kpp_const_fields%dto/kpp_const_fields%spd
  ocnT_time=(ndays_upd_ocnT)*(FLOOR(kpp_const_fields%time)*NINT(kpp_const_fields%spd)/&
       (kpp_const_fields%ndtupdocnT*NINT(kpp_const_fields%dto)))+&
       (0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdocnT)
  WRITE(nuout,*) ocnT_time,last_timein
     
  IF (ocnT_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_OCNT) THEN 
        DO WHILE (ocnT_time .gt. last_timein)
           ocnT_time=ocnT_time-kpp_const_fields%ocnT_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'Time for which to read the ocean temperatures exceeds the last time in the netCDF file &
             & and L_PERIODIC_OCNT has not been specified.  Attempting to read ocean temperatures will lead to &
             & an error, so aborting now ...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF

  write(nuout,*) 'Reading ocean temperature for time ',ocnT_time
  start(4)=NINT((ocnT_time-first_timein)*kpp_const_fields%spd/(kpp_const_fields%dto*kpp_const_fields%ndtupdocnT))+1
  
  write(nuout,*) 'Ocean temperature values are being read from position',start(4)
  status=NF_GET_VAR1_REAL(ocnT_ncid,time_varid,start(4),time_in)
  
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-ocnT_time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     write(nuerr,*) 'Cannot find time ',ocnT_time,' in ocean temperature climatology input file'
     write(nuerr,*) 'The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
  status=NF_GET_VARA_REAL(ocnT_ncid,ocnT_varid,start,count,ocnT_in)
  write(nuout,*) 'Ocean temperature climatology data have been read from position',start(4)
  
  ! Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
  ! into one long array with dimension NPTS.       
  DO ix=1,NX
     DO iy=1,NY
        ipoint=(iy-1)*nx+ix
        DO k=1,NZP1
           kpp_3d_fields%ocnT_clim(ipoint,k)=ocnT_in(ix,iy,k,1)
        ENDDO
     ENDDO
  ENDDO
  
  status=NF_CLOSE(ocnT_ncid)
  deallocate(ocnT_in)
  deallocate(longitudes)
  deallocate(latitudes)
  deallocate(z)
  
  RETURN
END SUBROUTINE MCKPP_READ_TEMPERATURES_3D

SUBROUTINE MCKPP_READ_TEMPERATURES_BOTTOM(kpp_3d_fields,kpp_const_fields)
  IMPLICIT NONE

  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER status,ncid
  REAL*4 var_in(NX,NY,1),time_in,first_timein, bottomclim_time,latitudes(NY_GLOBE),&
       longitudes(NX_GLOBE), last_timein, offset_temp
  INTEGER varid,time_varid,lat_varid,lon_varid,time_dimid,lat_dimid,lon_dimid,&
       nlat_file,nlon_file,ntime_file

  INTEGER count(3),start(3)
  INTEGER ix,iy,ipoint
  CHARACTER(LEN=30) tmp_name
  
  ! Set start and count to read a regional field.
  count(1)=NX
  count(2)=NY
  count(3)=1
  start(1)=1
  start(2)=1
  start(3)=1
  
  ! Open the netCDF file and find the latitude and longitude 
  ! boundaries in the input file.
  status=NF_OPEN(kpp_const_fields%bottom_file,0,ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  status=NF_INQ_VARID(ncid,'T',varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(ncid,'bottom temp climatology','latitude','longitude',&
       't',kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
  
  bottomclim_time=kpp_const_fields%time+0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdbottom
  IF (bottomclim_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_BOTTOM_TEMP) THEN 
        DO WHILE (bottomclim_time .gt. last_timein)
           bottomclim_time=bottomclim_time-kpp_const_fields%bottom_temp_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'Time for which to read bottom temperature exceeds the last time in the netCDF file and &
             & L_PERIODIC_BOTTOM_TEMP has not been specified.  Attempting to read bottom temperature will lead &
             & to an error, so aborting now...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF
 
  write(nuout,*) 'Reading climatological bottom temp for time ',bottomclim_time
  start(3)=NINT((bottomclim_time-first_timein)*kpp_const_fields%spd/&
       (kpp_const_fields%dto*kpp_const_fields%ndtupdbottom))+1
  write(nuout,*) 'Bottom temperatures are being read from position',start(3)
  write(nuout,*) 'start=',start,'count=',count
  
  status=NF_GET_VAR1_REAL(ncid,time_varid,start(3),time_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-bottomclim_time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     write(nuerr,*) 'Cannot find time',bottomclim_time,'in bottom temperature climatology file'
     write(nuerr,*) 'The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
  status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  WRITE(nuout,*) 'Bottom temperatures have been read from position',start(3)
  WRITE(nuout,*) 'First column of bottom temperatures: ',var_in(1,:,1)
  
!     KPP expects temperatures in CELSIUS.  If climatological bottom 
!     temperatures are in Kelvin, subtract 273.15.     
  offset_temp=0.
  ix=1
  iy=1
  DO WHILE (offset_temp.EQ.0.AND.ix.LE.NX)
     DO iy=1,NY
        IF (var_in(ix,iy,1) .gt. 200 .and. var_in(ix,iy,1) .lt. 400) offset_temp = 273.15    
     END DO
     ix=ix+1
  ENDDO
  
  DO ix=1,NX
     DO iy=1,NY
        ipoint=(iy-1)*NX+ix
        kpp_3d_fields%bottom_temp(ipoint) = var_in(ix,iy,1)-offset_temp        
     ENDDO
  ENDDO
  
  status=NF_CLOSE(ncid)
  
  RETURN
END SUBROUTINE MCKPP_READ_TEMPERATURES_BOTTOM
