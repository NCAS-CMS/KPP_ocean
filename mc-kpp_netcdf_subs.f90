SUBROUTINE MCKPP_NCDF_DEF_DIM (ncid,dimid,dimlen,varid,name,units,delta,long_name)
  
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
      
  INTEGER ncid,dimid,varid,dimlen
  REAL*4 delta
  CHARACTER*(*) name,units,long_name

  INTEGER status,unit_len,long_len

  unit_len=len(units)
  long_len=len(long_name)
  
  status=NF_DEF_DIM(ncid,name,dimlen,dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_DEF_VAR(ncid,name,NF_FLOAT,1,dimid,varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_PUT_ATT_TEXT (ncid,varid,'units',unit_len,units)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_PUT_ATT_REAL (ncid,varid,'spacing',NF_FLOAT,1,delta)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_PUT_ATT_TEXT (ncid,varid,'long_name',long_len,long_name)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  RETURN
END SUBROUTINE MCKPP_NCDF_DEF_DIM

SUBROUTINE MCKPP_NCDF_DEF_VAR (ncid,varid,ndims,dims,name,units,long_name)

  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  include 'netcdf.inc'
  
  REAL*4 valid_min,valid_max,miss_val,fill_val
  PARAMETER (valid_max=1.e19,valid_min=-valid_max,miss_val=1.e20,fill_val=miss_val)
  
  INTEGER ncid,varid,ndims
  INTEGER dims(ndims)
  CHARACTER*(*) name,units,long_name
  
  INTEGER status,unit_len,long_len
  
  unit_len=len(units)
  long_len=len(long_name)
  
  status=NF_DEF_VAR(ncid,name,nf_float,ndims,dims,varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_PUT_ATT_TEXT (ncid,varid,'long_name',long_len,long_name)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_PUT_ATT_TEXT (ncid,varid,'units',unit_len,units)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_PUT_ATT_REAL (ncid,varid,'valid_min',nf_float,1,valid_min)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_PUT_ATT_REAL (ncid,varid,'valid_max',nf_float,1,valid_max)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_PUT_ATT_REAL (ncid,varid,'missing_value',nf_float,1,miss_val)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_PUT_ATT_REAL (ncid,varid,'_FillValue',nf_float,1,fill_val)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  RETURN
END SUBROUTINE MCKPP_NCDF_DEF_VAR

SUBROUTINE MCKPP_HANDLE_ERR(status)
  
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  include 'netcdf.inc'
  
  integer status
  CHARACTER*80 err_message
  
  err_message=NF_STRERROR(status)
  write(nuerr,*) err_message
  
  CALL MCKPP_ABORT
END SUBROUTINE MCKPP_HANDLE_ERR

SUBROUTINE MCKPP_READ_PAR (kpp_3d_fields,ncid,vname,npars,nt,par_out)
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  
  CHARACTER(LEN=*) vname
  INTEGER ncid,npars,nt
  REAL par_out(npts,npars,nt)
  REAL*4,allocatable :: par_in(:,:,:,:),x_in(:),y_in(:)

  INTEGER start(4),count(4),ix,iy,ipt,ipar,ixx,iyy
  INTEGER status,dimid,varid
  
  allocate(par_in(nx,ny,npars,nt))
  allocate(x_in(NX_GLOBE))
  allocate(y_in(NY_GLOBE))
  
  status=NF_INQ_DIMID(ncid,'longitude',dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(ncid,'longitude',varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  WRITE(6,*) 'Routine read_par is reading longitude'
  status=NF_GET_VAR_REAL(ncid,varid,x_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ixx=1
  DO WHILE (abs(x_in(ixx)-kpp_3d_fields%dlon(1)) .GT. 1.e-3)
     ixx=ixx+1         
  ENDDO
  start(1)=ixx
  count(1)=nx
  
  status=NF_INQ_DIMID(ncid,'latitude',dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(ncid,'latitude',varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  WRITE(6,*) 'Routine read_par is reading latitude'
  status=NF_GET_VAR_REAL(ncid,varid,y_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  iyy=1
  DO WHILE (abs(y_in(iyy)-kpp_3d_fields%dlat(1)) .GT. 1.e-3)
     WRITE(6,*) y_in(iyy),kpp_3d_fields%dlat(1)
     iyy=iyy+1
  ENDDO
  start(2)=iyy
  count(2)=ny
  
  IF (vname .EQ. 'lsm') THEN
     DO iy=1,ny 
        DO ix=1,nx
           ipt=(iy-1)*nx+ix
           kpp_3d_fields%dlat(ipt)=y_in(iy+iyy-1)
           kpp_3d_fields%dlon(ipt)=x_in(ix+ixx-1) 
        ENDDO
     ENDDO
  ENDIF
  
  start(3)=1
  count(3)=npars
  start(4)=1
  count(4)=1
  
  status=NF_INQ_VARID(ncid,vname,varid)
  WRITE(6,*) 'Routine read_par is reading variable ',vname
  WRITE(6,*) 'start=',start,'count=',count
  status=NF_GET_VARA_REAL(ncid,varid,start,count,par_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  DO iy=1,ny
     DO ix=1,nx
        ipt=(iy-1)*nx+ix
        DO ipar=1,npars
           par_out(ipt,ipar,1)=par_in(ix,iy,ipar,1)
        ENDDO
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE MCKPP_READ_PAR

SUBROUTINE MCKPP_READ_IPAR (kpp_3d_fields,ncid,vname,npars,nt,par_out)
  
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
!Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  
  CHARACTER*20 vname
  INTEGER ncid,npars,nt
  INTEGER par_out(npts,npars,nt)
  INTEGER par_in(nx,ny,npars,nt)
  REAL*4 x_in(NX_GLOBE),y_in(NY_GLOBE)
  
  INTEGER start(4),count(4),ix,iy,ipt,ipar
  INTEGER status,dimid,varid
  
  status=NF_INQ_DIMID(ncid,'longitude',dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(ncid,'longitude',varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VAR_REAL(ncid,varid,x_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ix=1
  DO WHILE (abs(x_in(ix)-kpp_3d_fields%dlon(1)) .GT. 1.e-3)
     ix=ix+1
  ENDDO
  start(1)=ix
  count(1)=nx
  
  status=NF_INQ_DIMID(ncid,'latitude',dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(ncid,'latitude',varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VAR_REAL(ncid,varid,y_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  iy=1
  DO WHILE (abs(y_in(iy)-kpp_3d_fields%dlat(1)) .GT. 1.e-3)
     iy=iy+1
  ENDDO
  start(2)=iy
  count(2)=ny
  
  start(3)=1
  count(3)=npars
  start(4)=1
  count(4)=1
    
  status=NF_INQ_VARID(ncid,vname,varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VARA_INT(ncid,varid,start,count,par_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)         
  DO iy=1,ny
     DO ix=1,nx
        ipt=(iy-1)*nx+ix
        DO ipar=1,npars
           par_out(ipt,ipar,1)=par_in(ix,iy,ipar,1)
        ENDDO
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE MCKPP_READ_IPAR

SUBROUTINE MCKPP_DETERMINE_NETCDF_BOUNDARIES(ncid,file_description,latitude_name,longitude_name,&
     time_name,start_lon,start_lat,offset_lon,offset_lat,first_time,last_time,time_varid)

  IMPLICIT NONE
  
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <parameter.inc>
#include <netcdf.inc>

  INTEGER ncid,offset_lon,offset_lat,lon_dimid,lon_varid,lat_dimid,lat_varid,time_dimid,time_varid,ix,iy
  REAL start_lon,start_lat,first_time,last_time      
  CHARACTER(*) file_description,latitude_name,longitude_name,time_name
  CHARACTER(LEN=30) tmp_name
  
  INTEGER nlat_file,nlon_file,ntime_file,status
  REAL*4 longitudes(NX_GLOBE),latitudes(NY_GLOBE)
  
  status=NF_INQ_DIMID(ncid,latitude_name,lat_dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIM(ncid,lat_dimid,tmp_name,nlat_file)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)      
  status=NF_INQ_VARID(ncid,latitude_name,lat_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  WRITE(nuout,*) 'Inquired on '//latitude_name,lat_varid
  
  status=NF_INQ_DIMID(ncid,longitude_name,lon_dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIM(ncid,lon_dimid,tmp_name,nlon_file)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)      
  status=NF_INQ_VARID(ncid,longitude_name,lon_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  WRITE(nuout,*) 'Inquired on '//longitude_name,lon_varid
  
  status=NF_GET_VAR_REAL(ncid,lat_varid,latitudes(1:nlat_file))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  WRITE(nuout,*) 'Read ',nlat_file,' latitudes'
  status=NF_GET_VAR_REAL(ncid,lon_varid,longitudes(1:nlon_file))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  WRITE(nuout,*) 'Read ',nlon_file,' longitudes'
  
  ix=1
  DO WHILE (abs(longitudes(ix)-start_lon) .GT. 1.e-3)
     IF (ix .gt. nlon_file) THEN
        WRITE(nuout,*) 'Could not find starting longitude in the '//file_description//' file'
        CALL MCKPP_ABORT
     ENDIF
     ix=ix+1
  ENDDO
  offset_lon=ix
  
  iy=1
  DO WHILE (abs(latitudes(iy)-start_lat) .GT. 1.e-3)
     IF (iy .gt. nlat_file) THEN
        WRITE(nuout,*) 'Could not find starting latitude in the '//file_description//' file'
        CALL MCKPP_ABORT
     ENDIF
     iy=iy+1
  ENDDO
  offset_lat=iy
  
  status=NF_INQ_VARID(ncid,time_name,time_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VAR1_REAL(ncid,time_varid,1,first_time)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIMID(ncid,time_name,time_dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIM(ncid,time_dimid,tmp_name,ntime_file)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VAR1_REAL(ncid,time_varid,ntime_file,last_time)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  RETURN
END SUBROUTINE MCKPP_DETERMINE_NETCDF_BOUNDARIES
