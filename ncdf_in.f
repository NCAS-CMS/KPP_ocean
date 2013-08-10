      SUBROUTINE read_init (kpp_3d_fields,kpp_const_fields)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
      include 'vert_pgrid.com'
      include 'initialcon.com'
c      include 'location.com'
      include 'constants.com'
      
      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

      INTEGER status,ncid
      REAL*4, allocatable :: var_in(:,:,:),z_in(:),x_in(:),y_in(:)

      INTEGER varid,dimid
      INTEGER nz_in,nx_in,ny_in

      INTEGER ix,iy,ipt,count(3),start(3)
      INTEGER kin,k
      REAL deltaz,deltavar,offset_sst
      
      allocate(var_in(NX,NY,200))
      allocate(z_in(200))
      allocate(x_in(NX_GLOBE))
      allocate(y_in(NY_GLOBE))

      status=NF_OPEN(initdata_file,0,ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      status=NF_INQ_DIMID(ncid,'longitude',dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIMLEN (ncid,dimid,nx_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid,'longitude',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR_REAL(ncid,varid,x_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      ix=1
      DO WHILE (abs(x_in(ix)-kpp_3d_fields%dlon(1)) .GT. 1.e-3)
         ix=ix+1
         IF (ix .GE. nx_in) THEN
            write(nuerr,*) 'Error reading initial conditions'
            write(nuerr,*) 'Can''t find longitude ',
     +           kpp_3d_fields%dlon(1),' in range ',x_in(1),x_in(nx_in)
            CALL MIXED_ABORT
         ENDIF
      ENDDO
      start(1)=ix
      count(1)=nx

      status=NF_INQ_DIMID(ncid,'latitude',dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIMLEN (ncid,dimid,ny_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid,'latitude',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR_REAL(ncid,varid,y_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      iy=1
      DO WHILE (abs(y_in(iy)-kpp_3d_fields%dlat(1)) .GT. 1.e-3)
         iy=iy+1
         IF (iy .GE. ny_in) THEN
            write(nuerr,*) 'Error reading initial conditions'
            write(nuerr,*) 'Can''t find latitude ',
     +           kpp_3d_fields%dlat(1),' in range ',y_in(1),y_in(ny_in)
            CALL MIXED_ABORT
         ENDIF
      ENDDO
      start(2)=iy
      count(2)=ny

      status=NF_INQ_DIMID(ncid,'zvel',dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIMLEN (ncid,dimid,nz_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      start(3)=1
      count(3)=nz_in

      status=NF_INQ_VARID(ncid,'zvel',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR_REAL(ncid,varid,z_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid,'u',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      write(nuout,*) 'read_init read u'

      IF (L_INTERPINIT) THEN
         DO iy=1,ny
            DO ix=1,nx
               ipt=(iy-1)*nx+ix
               kin=1
               DO k=1,NZP1
                  IF (kpp_const_fields%zm(k) .GT. z_in(1)) THEN
                     kpp_3d_fields%U(ipt,k,1)=var_in(ix,iy,1)
                  ELSEIF (kpp_const_fields%zm(k) .LT. z_in(nz_in)) THEN
                     kpp_3d_fields%U(ipt,k,1)=var_in(ix,iy,nz_in)
                  ELSE               
                     DO WHILE (z_in(kin+1) .GT. kpp_const_fields%zm(k))
                        kin=kin+1
                     ENDDO
                     deltaz=z_in(kin)-z_in(kin+1)
                     deltavar=var_in(ix,iy,kin)-var_in(ix,iy,kin+1)
                     kpp_3d_fields%U(ipt,k,1)=var_in(ix,iy,kin)+
     &                    deltavar*(kpp_const_fields%zm(k)-z_in(kin))/
     +                    deltaz
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
c         U(ipt,:,1)=0.
      ELSE
         write(nuerr,*) 'You have to interpolate'
      ENDIF
      write(6,*) 'read_init interpolated u'
      status=NF_INQ_VARID(ncid,'v',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'read_init read v'

      IF (L_INTERPINIT) THEN
         DO iy=1,ny
            DO ix=1,nx
               ipt=(iy-1)*nx+ix
               kin=1
               DO k=1,NZP1
                  IF (kpp_const_fields%zm(k) .GT. z_in(1)) THEN
                     kpp_3d_fields%U(ipt,k,2)=var_in(ix,iy,1)
                  ELSEIF (kpp_const_fields%zm(k) .LT. z_in(nz_in)) THEN
                     kpp_3d_fields%U(ipt,k,2)=var_in(ix,iy,nz_in)
                  ELSE               
                     DO WHILE (z_in(kin+1) .GT. kpp_const_fields%zm(k))
                        kin=kin+1
                     ENDDO
                     deltaz=z_in(kin)-z_in(kin+1)
                     deltavar=var_in(ix,iy,kin)-var_in(ix,iy,kin+1)
                     kpp_3d_fields%U(ipt,k,2)=var_in(ix,iy,kin)+
     &                    deltavar*(kpp_const_fields%zm(k)-z_in(kin))/
     +                    deltaz
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
c         U(ipt,:,2)=0.
      ELSE
         write(nuerr,*) 'You have to interpolate'
      ENDIF
      write(6,*) 'read_init interpolated v'
c     Save initial currents in case they are needed to reinitalise
c     dodgy profiles (see resetting routines in steves_3d_ocn.f)
      kpp_3d_fields%U_init=kpp_3d_fields%U

      status=NF_INQ_DIMID(ncid,'ztemp',dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIMLEN (ncid,dimid,nz_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      start(3)=1
      count(3)=nz_in
      
      status=NF_INQ_VARID(ncid,'ztemp',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR_REAL(ncid,varid,z_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid,'temp',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'read_init read temp'
      
      IF (L_INTERPINIT) THEN
         DO iy=1,ny
            DO ix=1,nx
               ipt=(iy-1)*nx+ix
               kin=1
               DO k=1,NZP1
                  IF (kpp_const_fields%zm(k) .GT. z_in(1)) THEN
                     kpp_3d_fields%X(ipt,k,1)=var_in(ix,iy,1)
                  ELSEIF (kpp_const_fields%zm(k) .LT. z_in(nz_in)) THEN
                     kpp_3d_fields%X(ipt,k,1)=var_in(ix,iy,nz_in)
                  ELSE               
                     DO WHILE (z_in(kin+1) .GT. kpp_const_fields%zm(k))
                        kin=kin+1
                     ENDDO
                     deltaz=z_in(kin)-z_in(kin+1)
                     deltavar=var_in(ix,iy,kin)-var_in(ix,iy,kin+1)
                     kpp_3d_fields%X(ipt,k,1)=var_in(ix,iy,kin)+
     &                    deltavar*(kpp_const_fields%zm(k)-z_in(kin))/
     +                    deltaz
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ELSE
         write(nuerr,*) 'You have to interpolate'
      ENDIF
c
c     KPP requires temperatures in CELSIUS.  If initial conditions
c     are in Kelvin, subtract 273.15
c
      offset_sst = 0.
      DO ix=1,nx
         DO iy=1,ny
            IF (var_in(ix,iy,1) .gt. 200 .and. 
     &           var_in(ix,iy,1) .lt. 400)
     &           offset_sst = kpp_const_fields%TK0
         END DO 
      END DO
      kpp_3d_fields%X(:,:,1) = kpp_3d_fields%X(:,:,1) - offset_sst

      status=NF_INQ_DIMID(ncid,'zsal',dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIMLEN (ncid,dimid,nz_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      start(3)=1
      count(3)=nz_in

      status=NF_INQ_VARID(ncid,'zsal',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR_REAL(ncid,varid,z_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid,'sal',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'read_init read sal'
      
      IF (L_INTERPINIT) THEN
         DO iy=1,ny
            DO ix=1,nx
               ipt=(iy-1)*nx+ix               
               kin=1
               DO k=1,NZP1
                  IF (kpp_const_fields%zm(k) .GT. z_in(1)) THEN
                     kpp_3d_fields%X(ipt,k,2)=var_in(ix,iy,1)
                  ELSEIF (kpp_const_fields%zm(k) .LT. z_in(nz_in)) THEN
                     kpp_3d_fields%X(ipt,k,2)=var_in(ix,iy,nz_in)
                  ELSE               
                     DO WHILE (z_in(kin+1) .GT. kpp_const_fields%zm(k))
                        kin=kin+1
                     ENDDO
                     deltaz=z_in(kin)-z_in(kin+1)
                     deltavar=var_in(ix,iy,kin)-var_in(ix,iy,kin+1)
                     kpp_3d_fields%X(ipt,k,2)=var_in(ix,iy,kin)+
     &                    deltavar*(kpp_const_fields%zm(k)-z_in(kin))/
     +                    deltaz
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ELSE
         write(nuerr,*) 'You have to interpolate'
      ENDIF
      

      RETURN
      END

      SUBROUTINE init_flxdata(fname)
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <netcdf.inc>
#include <flx_in.com>
      
      CHARACTER*40 fname
      INTEGER status,index(3)
      
      index(1)=1
      index(2)=1
      index(3)=1

      status=NF_OPEN(fname,0,ncid_flx)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      status=NF_INQ_VARID(ncid_flx,'time',timein_id)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid_flx,'taux',varin_id(1))
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid_flx,'tauy',varin_id(2))
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid_flx,'swf',varin_id(3))
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid_flx,'lwf',varin_id(4))
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid_flx,'lhf',varin_id(5))
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid_flx,'shf',varin_id(6))
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid_flx,'precip',varin_id(7))
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      status=NF_GET_VAR1_REAL(ncid_flx,
     &     timein_id,index,first_timein)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

            
      RETURN
      END
      
      SUBROUTINE read_fluxes(taux,tauy,swf,lwf,lhf,shf,rain,snow,
     +     kpp_3d_fields,kpp_const_fields)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
#include <flx_in.com>

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      REAL taux(NPTS),tauy(NPTS),swf(NPTS),lwf(NPTS),
     $     lhf(NPTS),shf(NPTS),rain(NPTS),snow(NPTS)
c      REAL*4 var_in(NX,NY),time_in,time
      REAL*4 time_in,time
      
      INTEGER dimid,varid
      INTEGER ipt,ix,iy,nx_in,ny_in

c      REAL*4 x_in(NX_GLOBE),y_in(NY_GLOBE)
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

      status=NF_INQ_DIMID(ncid_flx,'longitude',dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIMLEN (ncid_flx,dimid,nx_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid_flx,'longitude',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR_REAL(ncid_flx,varid,x_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      ix=1
      DO WHILE (abs(x_in(ix)-kpp_3d_fields%dlon(1)) .GT. 1.e-3)
         ix=ix+1
         IF (ix .GE. nx_in) THEN
            write(nuerr,*) 'Error reading fluxes'
            write(nuerr,*) 'Can''t find longitude ',
     +           kpp_3d_fields%dlon(1),' in range ',x_in(1),x_in(nx_in)
            CALL MIXED_ABORT
         ENDIF
      ENDDO
      start(1)=ix
      count(1)=nx

      status=NF_INQ_DIMID(ncid_flx,'latitude',dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIMLEN (ncid_flx,dimid,ny_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid_flx,'latitude',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR_REAL(ncid_flx,varid,y_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      iy=1
      DO WHILE (abs(y_in(iy)-kpp_3d_fields%dlat(1)) .GT. 1.e-3)
         iy=iy+1
         IF (iy .GE. ny_in) THEN
            write(nuerr,*) 'Error reading fluxes'
            write(nuerr,*) 'Can''t find latitude ',
     +           kpp_3d_fields%dlat(1),' in range ',y_in(1),y_in(ny_in)
            CALL MIXED_ABORT
         ENDIF
      ENDDO
      start(2)=iy
      count(2)=ny


      time=kpp_const_fields%time+0.5*kpp_const_fields%dtsec/
     +     kpp_const_fields%spd
      write(nuout,*) ' Reading fluxes for time ',time,
     +     kpp_const_fields%time,kpp_const_fields%dtsec,
     +     kpp_const_fields%spd
c     IF (ndtocn .EQ. 1) THEN 
      start(3)=MAX(NINT((time-first_timein)*kpp_const_fields%spd/
     +     kpp_const_fields%dtsec)+1,1)
      WRITE(6,*) 'Reading time from time point ',start(3)
      WRITE(6,*) 'first_timein=',first_timein,'time=',time
      status=NF_GET_VAR1_REAL(ncid_flx,timein_id,start(3),time_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
#ifndef OASIS3      
      IF (abs(time_in-time) .GT. 0.01*kpp_const_fields%dtsec/
     +     kpp_const_fields%spd) THEN
         write(nuerr,*) 'Cannot find time,',time,'in fluxes file'
         write(nuerr,*) 'The closest I came was',time_in
         CALL MIXED_ABORT
      ENDIF
#endif
      WRITE(6,*) 'Reading fluxes from time point ',start(3)
      status=NF_GET_VARA_REAL(ncid_flx,varin_id(1),start,count
     $     ,var_in)
      WRITE(6,*) 'Read taux'
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      DO iy=1,ny
         DO ix=1,nx
            ipt=(iy-1)*nx+ix            
            taux(ipt)=var_in(ix,iy)
         ENDDO
      ENDDO      
      status=NF_GET_VARA_REAL(ncid_flx,varin_id(2),start,count
     $     ,var_in)
      WRITE(6,*) 'Read tauy'
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      DO iy=1,ny
         DO ix=1,nx
            ipt=(iy-1)*nx+ix
            tauy(ipt)=var_in(ix,iy)
         ENDDO
      ENDDO
      status=NF_GET_VARA_REAL(ncid_flx,varin_id(3),start,count
     $     ,var_in)
      WRITE(6,*) 'Read swf'
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      DO iy=1,ny
         DO ix=1,nx
            ipt=(iy-1)*nx+ix
            swf(ipt)=var_in(ix,iy)
         ENDDO
      ENDDO
      status=NF_GET_VARA_REAL(ncid_flx,varin_id(4),start,count
     $     ,var_in)
      WRITE(6,*) 'Read lwf'
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      DO iy=1,ny
         DO ix=1,nx
            ipt=(iy-1)*nx+ix
            lwf(ipt)=var_in(ix,iy)
         ENDDO
      ENDDO
      status=NF_GET_VARA_REAL(ncid_flx,varin_id(5),start,count
     $     ,var_in)
      WRITE(6,*) 'Read lhf'
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      DO iy=1,ny
         DO ix=1,nx
            ipt=(iy-1)*nx+ix
            lhf(ipt)=var_in(ix,iy)
         ENDDO
      ENDDO
      status=NF_GET_VARA_REAL(ncid_flx,varin_id(6),start,count
     $     ,var_in)
      WRITE(6,*) 'Read shf'
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)      
      DO iy=1,ny
         DO ix=1,nx
            ipt=(iy-1)*nx+ix
            shf(ipt)=var_in(ix,iy)
         ENDDO
      ENDDO
      status=NF_GET_VARA_REAL(ncid_flx,varin_id(7),start,count
     $     ,var_in)
      WRITE(6,*) 'Read rain'
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      DO iy=1,ny
         DO ix=1,nx
            ipt=(iy-1)*nx+ix
            rain(ipt)=var_in(ix,iy)
         ENDDO
      ENDDO
      DO iy=1,ny
         DO ix=1,nx
            ipt=(iy-1)*nx+ix
            snow(ipt)=0.0
         ENDDO
      ENDDO
      WRITE(6,*) 'Set snow to zero'
c     ELSE 
c     write(nuerr,*) 'You need some more code'
c     CALL MIXED_ABORT
c     ENDIF
c      WRITE(6,*) 'Setting time equal to dummy_time'
c      time=dummy_time
      
      WRITE(6,*) 'Finished reading fluxes'

      RETURN
      END

      SUBROUTINE init_parasfile
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
#include <parameter.inc>
      include 'proc_pars.com'
      
      INTEGER status      

      status=NF_OPEN(paras_file,0,ncid_paras)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      RETURN
      END

      SUBROUTINE close_parasfile

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
#include <parameter.inc>
      include 'proc_pars.com'
      
      INTEGER status

      status=NF_CLOSE(ncid_paras)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      RETURN
      END

      SUBROUTINE init_advectfile
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
#include <parameter.inc>
      include 'ocn_advec.com'
      
      INTEGER status      

      status=NF_OPEN(advect_file,0,ncid_advec)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      RETURN
      END

      SUBROUTINE close_advectfile

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
#include <parameter.inc>
      include 'ocn_advec.com'
      
      INTEGER status

      status=NF_CLOSE(ncid_advec)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      RETURN
      END

      SUBROUTINE init_landseafile
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
#include <parameter.inc>
#include <landsea.com>

      INTEGER status      
      
      WRITE(nuout,*) 'In init_landseafile'
      status=NF_OPEN(landsea_file,0,ncid_landsea)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      RETURN
      END

      SUBROUTINE close_landseafile

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
#include <parameter.inc>
#include <landsea.com>

      
      INTEGER status

      status=NF_CLOSE(ncid_landsea)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      RETURN
      END

      SUBROUTINE init_cplwghtfile
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER(nuout=6,nuerr=0)
      
      include 'netcdf.inc'
#include <parameter.inc>
      include 'couple.com'

      INTEGER status
     
      status=NF_OPEN(cplwght_file,0,ncid_cplwght)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      RETURN
      END
      
      SUBROUTINE close_cplwghtfile
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER(nuout=6,nuerr=0)

      include 'netcdf.inc'
#include <parameter.inc>
      include 'couple.com'

      INTEGER status

      status=NF_CLOSE(ncid_cplwght)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      RETURN
      END

      SUBROUTINE READ_PAR (kpp_3d_fields,
     &     ncid,vname,npars,nt,par_out)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
c      include 'location.com'

      TYPE(kpp_3d_type) :: kpp_3d_fields

      CHARACTER(LEN=*) vname
      INTEGER ncid,npars,nt
      REAL par_out(npts,npars,nt)
c      REAL*4 par_in(nx,ny,npars,nt)
c      REAL*4 x_in(NX_GLOBE),y_in(NY_GLOBE)
      REAL*4,allocatable :: par_in(:,:,:,:),x_in(:),y_in(:)

      INTEGER start(4),count(4),ix,iy,ipt,ipar,ixx,iyy
      INTEGER status,dimid,varid
      
      allocate(par_in(nx,ny,npars,nt))
      allocate(x_in(NX_GLOBE))
      allocate(y_in(NY_GLOBE))

      status=NF_INQ_DIMID(ncid,'longitude',dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid,'longitude',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(6,*) 'Routine read_par is reading longitude'
      status=NF_GET_VAR_REAL(ncid,varid,x_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      ixx=1
      DO WHILE (abs(x_in(ixx)-kpp_3d_fields%dlon(1)) .GT. 1.e-3)
         ixx=ixx+1         
      ENDDO
      start(1)=ixx
      count(1)=nx

      status=NF_INQ_DIMID(ncid,'latitude',dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid,'latitude',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(6,*) 'Routine read_par is reading latitude'
      status=NF_GET_VAR_REAL(ncid,varid,y_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
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
C               WRITE(nuout,*) 'ipt=',ipt,'dlat=',dlat(ipt)
C               WRITE(nuout,*) 'dlon=',dlon(ipt),'iyy=',iyy,'ixx=',ixx
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
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      DO iy=1,ny
         DO ix=1,nx
            ipt=(iy-1)*nx+ix
            DO ipar=1,npars
               par_out(ipt,ipar,1)=par_in(ix,iy,ipar,1)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE READ_IPAR (kpp_3d_fields,
     &     ncid,vname,npars,nt,par_out)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
!Automatically includes parameter.inc!
#include <kpp_3d_type.com>
c      include 'location.com'

      TYPE(kpp_3d_type) :: kpp_3d_fields

      CHARACTER*20 vname
      INTEGER ncid,npars,nt
      INTEGER par_out(npts,npars,nt)
      INTEGER par_in(nx,ny,npars,nt)
      REAL*4 x_in(NX_GLOBE),y_in(NY_GLOBE)

      INTEGER start(4),count(4),ix,iy,ipt,ipar
      INTEGER status,dimid,varid
      
      status=NF_INQ_DIMID(ncid,'longitude',dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid,'longitude',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR_REAL(ncid,varid,x_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      ix=1
      DO WHILE (abs(x_in(ix)-kpp_3d_fields%dlon(1)) .GT. 1.e-3)
         ix=ix+1
      ENDDO
      start(1)=ix
      count(1)=nx

      status=NF_INQ_DIMID(ncid,'latitude',dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid,'latitude',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR_REAL(ncid,varid,y_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
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
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VARA_INT(ncid,varid,start,count,par_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)         
      DO iy=1,ny
         DO ix=1,nx
            ipt=(iy-1)*nx+ix
            DO ipar=1,npars
               par_out(ipt,ipar,1)=par_in(ix,iy,ipar,1)
            ENDDO
         ENDDO
      ENDDO

      
      RETURN
      END

      SUBROUTINE read_fcorr(kpp_3d_fields,kpp_const_fields)

      IMPLICIT NONE
      INTEGER nuout,nuerr,start(3),count(3)
      INTEGER ix,iy,ipoint,fcorr_varid,status,lat_varid,lon_varid,
     +     time_varid,lat_dimid,lon_dimid,time_dimid,fcorr_ncid,k,
     +     nlat_file,nlon_file,ntime_file
      PARAMETER (nuout=6,nuerr=0)
      
      include 'netcdf.inc'
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
c      include 'location.com'
      include 'fcorr_in.com'
      include 'times.com'
      include 'timocn.com'
      include 'constants.com'

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      REAL*4 ixx,jyy,fcorr_twod_in(NX,NY,1),latitudes(NY),
     +     longitudes(NX),z(NZP1),first_timein,time_in,fcorr_time,
     +     ndays_upd_fcorr,last_timein
      CHARACTER(LEN=30) tmp_name

c     Read in a NetCDF file containing a time-varying flux correction
c     at the surface only.  Frequency of read is controlled by ndtupdfcorr
c     in the namelist
c
c     NPK 29/06/08

      status = NF_OPEN(fcorrin_file,0,fcorr_ncid)
      IF (status.NE.0) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Opened flux-correction input file'
      
      count(1)=NX
      count(2)=NY      
      count(3)=1
      start(1)=1
      start(2)=1
      start(3)=1
 
      CALL determine_netcdf_boundaries(fcorr_ncid,'flux correction',
     &     'latitude','longitude','t',kpp_3d_fields%dlon(1),
     +     kpp_3d_fields%dlat(1),start(1),
     &     start(2),first_timein,last_timein,time_varid)
 
      status=NF_INQ_VARID(fcorr_ncid,'fcorr',fcorr_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
       
      ndays_upd_fcorr = ndtupdfcorr*kpp_const_fields%dto/
     +     kpp_const_fields%spd
      fcorr_time=(ndays_upd_fcorr)*
     &     (FLOOR(kpp_const_fields%time)*NINT(kpp_const_fields%spd)/
     +     (ndtupdfcorr*NINT(kpp_const_fields%dto)))+
     &     (0.5*kpp_const_fields%dto/kpp_const_fields%spd*ndtupdfcorr)

      IF (fcorr_time .gt. last_timein) THEN
         IF (L_PERIODIC_FCORR) THEN 
            DO WHILE (fcorr_time .gt. last_timein)
               fcorr_time=fcorr_time-fcorr_period
            ENDDO
         ELSE
            WRITE(nuout,*) 'Time for which to read the flux & 
     & corrections exceeds the last time in the netCDF file &
     & and L_PERIODIC_FCORR has not been specified.  &
     & Attempting to read flux corrections will lead to an error, so &
     & aborting now ...'
            CALL MIXED_ABORT
         ENDIF
      ENDIF               

      write(nuout,*) 'Reading flux correction for time ',fcorr_time
      start(3)=NINT((fcorr_time-first_timein)*kpp_const_fields%spd/
     +     (kpp_const_fields%dto*ndtupdfcorr))+1
      write(nuout,*) 'Flux corrections are being read from position',
     &     start(3)
      status=NF_GET_VAR1_REAL(fcorr_ncid,time_varid,start(3),time_in)
      
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (abs(time_in-fcorr_time) .GT. 0.01*kpp_const_fields%dtsec/
     +     kpp_const_fields%spd) THEN
         write(nuerr,*) 'Cannot find time',fcorr_time,
     &        'in flux-correction input file'
         write(nuerr,*) 'The closest I came was',time_in
         CALL MIXED_ABORT
      ENDIF
      status=NF_GET_VARA_REAL(fcorr_ncid,fcorr_varid,start,count
     &     ,fcorr_twod_in)
      write(nuout,*) 'Flux corrections have been read from position',
     &     start(3)

c
c     Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
c     into one long array with dimension NPTS.
c         
      DO ix=1,NX
         DO iy=1,NY
            ipoint=(iy-1)*nx+ix
            kpp_3d_fields%fcorr_twod(ipoint)=fcorr_twod_in(ix,iy,1)
         ENDDO
      ENDDO

      status=NF_CLOSE(fcorr_ncid)

      END

      SUBROUTINE read_fcorrwithz(kpp_3d_fields,kpp_const_fields)
      
      IMPLICIT NONE
      INTEGER nuout,nuerr,start(4),count(4)
      INTEGER ix,iy,iz,ipoint,fcorr_varid,status,lat_varid,lon_varid,
     +     z_varid,z_dimid,time_varid,fcorr_ncid,k,lat_dimid,lon_dimid,
     +     time_dimid,nlon_file,nlat_file,ntime_file,nz_file

      PARAMETER (nuout=6,nuerr=0)
      include 'netcdf.inc'
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
c      include 'location.com'
      include 'fcorr_in.com'
      include 'times.com'
      include 'timocn.com'
      include 'constants.com'

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      REAL*4 ixx,jyy,first_timein,time_in,
     +     fcorr_time,ndays_upd_fcorr,last_timein
      CHARACTER(LEN=30) tmp_name
      REAL*4, allocatable :: fcorr_in(:,:,:,:),
     +     longitudes(:),latitudes(:),z(:)
c      COMMON /save_fcorr_withz/ fcorr_withz, Tinc_fcorr
c
c     Read in a NetCDF file containing a 
c     time-varying flux correction at every model vertical level.
c     Frequency of read is controlled by ndtupdfcorr in the namelist
c
c     NPK 12/02/08
c      
      allocate(fcorr_in(NX,NY,NZP1,1))
      allocate(longitudes(NX_GLOBE))
      allocate(latitudes(NY_GLOBE))
      allocate(z(NZP1))
      status=NF_OPEN(fcorrin_file,0,fcorr_ncid)
      IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Opened flux-correction input file ',fcorr_ncid
      
      count=(/NX,NY,NZP1,1/)
      start=(/1,1,1,1/)

      status=NF_INQ_VARID(fcorr_ncid,'z',z_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)  
      status=NF_INQ_DIMID(fcorr_ncid,'z',z_dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIM(fcorr_ncid,z_dimid,tmp_name,nz_file)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (NZP1.ne.nz_file) THEN
         WRITE(nuout,*) 'Input file for flux corrections does ',
     &        'not have the correct number of vertical levels. ',
     &        'It should have ',NZP1,' but instead has ',nz_file
         CALL MIXED_ABORT
      ELSE
         status=NF_GET_VAR_REAL(fcorr_ncid,z_varid,z)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         WRITE(nuout,*) 'Read in depths from the flux-correction ',
     &        'input file'
      ENDIF

      CALL determine_netcdf_boundaries(fcorr_ncid,'flux correction',
     &     'latitude','longitude','t',kpp_3d_fields%dlon(1),
     +     kpp_3d_fields%dlat(1),start(1),
     &     start(2),first_timein,last_timein,time_varid)

      status=NF_INQ_VARID(fcorr_ncid,'fcorr',fcorr_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      ndays_upd_fcorr = ndtupdfcorr*kpp_const_fields%dto/
     +     kpp_const_fields%spd
      WRITE(nuout,*) ndays_upd_fcorr,FLOOR(kpp_const_fields%time)*
     +     NINT(kpp_const_fields%spd),
     &     ndtupdfcorr*NINT(kpp_const_fields%dto),
     +     0.5*kpp_const_fields%dto/kpp_const_fields%spd*ndtupdfcorr
      fcorr_time=(ndays_upd_fcorr)*
     &     FLOOR(kpp_const_fields%time)*NINT(kpp_const_fields%spd)/
     +     FLOAT(ndtupdfcorr*NINT(kpp_const_fields%dto))+
     &     (0.5*kpp_const_fields%dto/kpp_const_fields%spd*ndtupdfcorr)
      WRITE(nuout,*) fcorr_time,last_timein
      
      IF (fcorr_time .gt. last_timein) THEN
         IF (L_PERIODIC_FCORR) THEN 
            DO WHILE (fcorr_time .gt. last_timein)
               fcorr_time=fcorr_time-fcorr_period
            ENDDO
         ELSE
            WRITE(nuout,*) 'Time for which to read the flux
     & corrections exceeds the last time in the netCDF file
     & and L_PERIODIC_FCORR has not been specified.
     & Attempting to read flux corrections will lead to an error, so
     & aborting now ...'
            CALL MIXED_ABORT
         ENDIF
      ENDIF               

      write(nuout,*) 'Reading flux correction for time ',fcorr_time
      start(4)=NINT((fcorr_time-first_timein)*kpp_const_fields%spd/
     +     (kpp_const_fields%dto*ndtupdfcorr))+1
      write(nuout,*) 'Flux corrections are being read from position',
     &     start(4)
      status=NF_GET_VAR1_REAL(fcorr_ncid,time_varid,start(4),time_in)
      
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (abs(time_in-fcorr_time) .GT. 0.01*kpp_const_fields%dtsec/
     +     kpp_const_fields%spd) THEN
         write(nuerr,*) 'Cannot find time',fcorr_time,
     &        'in flux-correction input file'
         write(nuerr,*) 'The closest I came was',time_in
         CALL MIXED_ABORT
      ENDIF
      status=NF_GET_VARA_REAL(fcorr_ncid,fcorr_varid,start,count
     &     ,fcorr_in)
      write(nuout,*) 'Flux corrections have been read from position',
     &     start(4)      
c
c     Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
c     into one long array with dimension NPTS.
c         
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
      END

!sfcorr added LH 24/05/2013

      SUBROUTINE read_sfcorr(kpp_3d_fields,kpp_const_fields)

      IMPLICIT NONE
      INTEGER nuout,nuerr,start(3),count(3)
      INTEGER ix,iy,ipoint,sfcorr_varid,status,lat_varid,lon_varid,
     +     time_varid,lat_dimid,lon_dimid,time_dimid,sfcorr_ncid,k,
     +     nlat_file,nlon_file,ntime_file
      PARAMETER (nuout=6,nuerr=0)
      
      include 'netcdf.inc'
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
c      include 'location.com'
      include 'sfcorr_in.com'
      include 'times.com'
      include 'timocn.com'
      include 'constants.com'

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      REAL*4 ixx,jyy,sfcorr_twod_in(NX,NY,1),latitudes(NY),
     +     longitudes(NX),z(NZP1),first_timein,time_in,sfcorr_time,
     +     ndays_upd_sfcorr,last_timein
c      COMMON /save_sfcorr_twod/ sfcorr_twod
      CHARACTER(LEN=30) tmp_name

c     Read in a NetCDF file containing a time-varying flux correction
c     at the surface only.  Frequency of read is controlled by ndtupdsfcorr
c     in the namelist
c
c     NPK 29/06/08

      status = NF_OPEN(sfcorrin_file,0,sfcorr_ncid)
      IF (status.NE.0) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Opened flux-correction input file'
      
      count(1)=NX
      count(2)=NY      
      count(3)=1
      start(1)=1
      start(2)=1
      start(3)=1
 
      CALL determine_netcdf_boundaries(sfcorr_ncid,'flux correction',
     &     'latitude','longitude','t',kpp_3d_fields%dlon(1),
     +     kpp_3d_fields%dlat(1),start(1),
     &     start(2),first_timein,last_timein,time_varid)
 
      status=NF_INQ_VARID(sfcorr_ncid,'sfcorr',sfcorr_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
       
      ndays_upd_sfcorr = ndtupdsfcorr*kpp_const_fields%dto/
     +     kpp_const_fields%spd
      sfcorr_time=(ndays_upd_sfcorr)*
     &     (FLOOR(kpp_const_fields%time)*NINT(kpp_const_fields%spd)/
     +     (ndtupdsfcorr*NINT(kpp_const_fields%dto)))+
     &     (0.5*kpp_const_fields%dto/kpp_const_fields%spd*ndtupdsfcorr)

      IF (sfcorr_time .gt. last_timein) THEN
         IF (L_PERIODIC_SFCORR) THEN 
            DO WHILE (sfcorr_time .gt. last_timein)
               sfcorr_time=sfcorr_time-sfcorr_period
            ENDDO
         ELSE
            WRITE(nuout,*) 'Time for which to read the flux & 
     & corrections exceeds the last time in the netCDF file &
     & and L_PERIODIC_SFCORR has not been specified.  &
     & Attempting to read flux corrections will lead to an error, so &
     & aborting now ...'
            CALL MIXED_ABORT
         ENDIF
      ENDIF               

      write(nuout,*) 'Reading flux correction for time ',sfcorr_time
      start(3)=NINT((sfcorr_time-first_timein)*kpp_const_fields%spd/
     +     (kpp_const_fields%dto*ndtupdsfcorr))+1
      write(nuout,*) 'Flux corrections are being read from position',
     &     start(3)
      status=NF_GET_VAR1_REAL(sfcorr_ncid,time_varid,start(3),time_in)
      
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (abs(time_in-sfcorr_time) .GT. 0.01*kpp_const_fields%dtsec/
     +     kpp_const_fields%spd) THEN
         write(nuerr,*) 'Cannot find time',sfcorr_time,
     &        'in flux-correction input file'
         write(nuerr,*) 'The closest I came was',time_in
         CALL MIXED_ABORT
      ENDIF
      status=NF_GET_VARA_REAL(sfcorr_ncid,sfcorr_varid,start,count
     &     ,sfcorr_twod_in)
      write(nuout,*) 'Flux corrections have been read from position',
     &     start(3)

c
c     Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
c     into one long array with dimension NPTS.
c         
      DO ix=1,NX
         DO iy=1,NY
            ipoint=(iy-1)*nx+ix
            kpp_3d_fields%sfcorr_twod(ipoint)=sfcorr_twod_in(ix,iy,1)
         ENDDO
      ENDDO

      status=NF_CLOSE(sfcorr_ncid)

      END

!sfcorr_withz added LH 24/05/2013

      SUBROUTINE read_sfcorrwithz(kpp_3d_fields,kpp_const_fields)
      
      IMPLICIT NONE
      INTEGER nuout,nuerr,start(4),count(4)
      INTEGER ix,iy,iz,ipoint,sfcorr_varid,status,lat_varid,lon_varid,
     +     z_varid,z_dimid,time_varid,sfcorr_ncid,k,lat_dimid,lon_dimid,
     +     time_dimid,nlon_file,nlat_file,ntime_file,nz_file

      PARAMETER (nuout=6,nuerr=0)
      include 'netcdf.inc'
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
c      include 'location.com'
      include 'sfcorr_in.com'
      include 'times.com'
      include 'timocn.com'
      include 'constants.com'

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      REAL*4 ixx,jyy,first_timein,time_in,
     +     sfcorr_time,ndays_upd_sfcorr,last_timein
      REAL*4,allocatable :: sfcorr_in(:,:,:,:),latitudes(:),
     +     longitudes(:),z(:)
      CHARACTER(LEN=30) tmp_name
c      COMMON /save_sfcorr_withz/ sfcorr_withz
c
c     Read in a NetCDF file containing a 
c     time-varying flux correction at every model vertical level.
c     Frequency of read is controlled by ndtupdsfcorr in the namelist
c
c     NPK 12/02/08
c      
      allocate(sfcorr_in(NX,NY,NZP1,1))
      allocate(latitudes(NX_GLOBE))
      allocate(longitudes(NY_GLOBE))
      allocate(z(NZP1))
      status=NF_OPEN(sfcorrin_file,0,sfcorr_ncid)
      IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Opened salinity-correction input file ',
     +     sfcorr_ncid
      
      count=(/NX,NY,NZP1,1/)
      start=(/1,1,1,1/)

      status=NF_INQ_VARID(sfcorr_ncid,'z',z_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)  
      status=NF_INQ_DIMID(sfcorr_ncid,'z',z_dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIM(sfcorr_ncid,z_dimid,tmp_name,nz_file)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (NZP1.ne.nz_file) THEN
         WRITE(nuout,*) 'Input file for salinity corrections does ',
     &        'not have the correct number of vertical levels. ',
     &        'It should have ',NZP1,' but instead has ',nz_file
         CALL MIXED_ABORT
      ELSE
         status=NF_GET_VAR_REAL(sfcorr_ncid,z_varid,z)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         WRITE(nuout,*) 'Read in depths from the salinity-correction ',
     &        'input file'
      ENDIF

      CALL determine_netcdf_boundaries(sfcorr_ncid,'salinity correction'
     &     ,'latitude','longitude','t',kpp_3d_fields%dlon(1),
     +     kpp_3d_fields%dlat(1),start(1),
     &     start(2),first_timein,last_timein,time_varid)

      status=NF_INQ_VARID(sfcorr_ncid,'sfcorr',sfcorr_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      ndays_upd_sfcorr = ndtupdsfcorr*kpp_const_fields%dto/
     +     kpp_const_fields%spd
      WRITE(nuout,*) ndays_upd_sfcorr,FLOOR(kpp_const_fields%time)*
     +     NINT(kpp_const_fields%spd),
     &     ndtupdsfcorr*NINT(kpp_const_fields%dto),
     +     0.5*kpp_const_fields%dto/kpp_const_fields%spd*ndtupdsfcorr
      sfcorr_time=(ndays_upd_sfcorr)*
     &     FLOOR(kpp_const_fields%time)*NINT(kpp_const_fields%spd)/
     +     FLOAT(ndtupdsfcorr*NINT(kpp_const_fields%dto))+
     &     (0.5*kpp_const_fields%dto/kpp_const_fields%spd*ndtupdsfcorr)
      WRITE(nuout,*) sfcorr_time,last_timein
      
      IF (sfcorr_time .gt. last_timein) THEN
         IF (L_PERIODIC_SFCORR) THEN 
            DO WHILE (sfcorr_time .gt. last_timein)
               sfcorr_time=sfcorr_time-sfcorr_period
            ENDDO
         ELSE
            WRITE(nuout,*) 'Time for which to read the flux
     & corrections exceeds the last time in the netCDF file
     & and L_PERIODIC_SFCORR has not been specified.
     & Attempting to read salinity corrections will lead to an error, so
     & aborting now ...'
            CALL MIXED_ABORT
         ENDIF
      ENDIF               

      write(nuout,*) 'Reading salinity correction for time ',sfcorr_time
      start(4)=NINT((sfcorr_time-first_timein)*kpp_const_fields%spd/
     +     (kpp_const_fields%dto*ndtupdsfcorr))+1
      write(nuout,*) 'salinity corrections are being read from position'
     &     ,start(4)
      status=NF_GET_VAR1_REAL(sfcorr_ncid,time_varid,start(4),time_in)
      
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (abs(time_in-sfcorr_time) .GT. 0.01*kpp_const_fields%dtsec/
     +     kpp_const_fields%spd) THEN
         write(nuerr,*) 'Cannot find time',sfcorr_time,
     &        'in flux-correction input file'
         write(nuerr,*) 'The closest I came was',time_in
         CALL MIXED_ABORT
      ENDIF
      status=NF_GET_VARA_REAL(sfcorr_ncid,sfcorr_varid,start,count
     &     ,sfcorr_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      write(nuout,*) 'salinity corrections have been read from position'
     &     ,start(4)
      
c
c     Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
c     into one long array with dimension NPTS.
c         
      DO ix=1,NX
         DO iy=1,NY
            ipoint=(iy-1)*nx+ix
            DO k=1,NZP1
               kpp_3d_fields%sfcorr_withz(ipoint,k)=sfcorr_in(ix,iy,k,1)
            ENDDO
         ENDDO
      ENDDO

      status=NF_CLOSE(sfcorr_ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)      
      deallocate(sfcorr_in)
      deallocate(longitudes)
      deallocate(latitudes)
      deallocate(z)
      
      RETURN
      END

#ifdef COUPLE
#ifdef CFS
      SUBROUTINE read_landsea_global
      IMPLICIT NONE
      
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
#include <parameter.inc>
#include <landsea.com>
      include 'constants.com'
c      include 'location.com'
      
      INTEGER status,varid,latid,lonid
      REAL*4 var_in(NX_GLOBE,NY_GLOBE),latitudes_in(NY_GLOBE),
     +     longitudes_in(NX_GLOBE)
      INTEGER count(3),start(3),ix,iy
      
      status=NF_OPEN(landsea_file,0,ncid_landsea)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)      
      
      count(1)=NX_GLOBE
      count(2)=NY_GLOBE
      count(3)=1

      start(1)=1
      start(2)=1
      start(3)=1
      
      status=NF_INQ_VARID(ncid_landsea,'latitude',latid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid_landsea,'longitude',lonid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR_REAL(ncid_landsea,latid,latitudes_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR_REAL(ncid_landsea,lonid,longitudes_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      global_longitudes=longitudes_in
      global_latitudes=latitudes_in

      status=NF_INQ_VARID(ncid_landsea,'lsm',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VARA_REAL(ncid_landsea,varid,start,count,var_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      global_lsm=var_in

      status=NF_CLOSE(ncid_landsea)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      RETURN
      END
#endif
#endif

      SUBROUTINE read_sstin(kpp_3d_fields,kpp_const_fields)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
      include 'constants.com'
      include 'couple.com'
      include 'times.com'
      include 'timocn.com'
      include 'sstclim.com'
c      include 'location.com'
#include <currclim.com>

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      INTEGER sst_nx,sst_ny
#ifdef COUPLE
      PARAMETER(sst_nx=NX_GLOBE,sst_ny=NY_GLOBE)
#else
      PARAMETER(sst_nx=NX,sst_ny=NY)
#endif
      REAL sst_in(sst_nx,sst_ny,1),ice_in(sst_nx,sst_ny,1),
     &     snowdepth_in(sst_nx,sst_ny,1),icedepth_in(sst_nx,sst_ny,1)
      REAL usf_in(sst_nx,sst_ny),vsf_in(sst_nx,sst_ny),
     &     offset_sst
      INTEGER status,ncid
      REAL*4 var_in(sst_nx,sst_ny,1),time_in,
     &     first_timein,last_timein,sstclim_time,longitudes(NX_GLOBE),
     &     latitudes(NY_GLOBE)
      INTEGER varid, time_varid,lat_varid,lon_varid,
     &     lon_dimid,lat_dimid,time_dimid
      INTEGER count(3),start(3)
      INTEGER ix,iy,nlat_file,nlon_file,ntime_file
      CHARACTER(LEN=30) tmp_name

      COMMON /save_sstin/ sst_in,ice_in,snowdepth_in,
     &     icedepth_in,usf_in,vsf_in

c      WRITE(nuout,*) 'Trying to open the sstin_file=',sstin_file
      status=NF_OPEN(sstin_file,0,ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
c      WRITE(nuout,*) 'Opened the sstin_file=',sstin_file

      count(1)=sst_nx
      count(2)=sst_ny
      count(3)=1
      start(1)=1
      start(2)=1
      start(3)=1

      status=NF_INQ_VARID(ncid,'sst',varid)
      IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)

#ifndef COUPLE
      CALL determine_netcdf_boundaries(ncid,'SST climatology',
     &     'latitude','longitude','time',kpp_3d_fields%dlon(1),
     +     kpp_3d_fields%dlat(1),start(1),
     &     start(2),first_timein,last_timein,time_varid)
#else
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid,'time',time_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIMID(ncid,'time',time_dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIM(ncid,time_dimid,tmp_name,ntime_file)
      status=NF_GET_VAR1_REAL(ncid,
     &     time_varid,start(3),first_timein)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR1_REAL(ncid,
     &     time_varid,ntime_file,last_timein)
#endif
      sstclim_time=kpp_const_fields%time+0.5*kpp_const_fields%dto/
     +     kpp_const_fields%spd*ndtupdsst
c      WRITE(6,*) kpp_const_fields%time,kpp_const_fields%dto,
c     +     kpp_const_fields%spd
      IF (sstclim_time .gt. last_timein) THEN
         IF (L_PERIODIC_CLIMSST) THEN 
            DO WHILE (sstclim_time .gt. last_timein)
               sstclim_time=sstclim_time-climsst_period
            ENDDO
         ELSE
            WRITE(nuout,*) 'Time for which to read SST exceeds
     &the last time in the netCDF file and L_PERIODIC_CLIMSST has
     &not been specified.  Attempting to read SST will lead to
     &an error, so aborting now ...'
            CALL MIXED_ABORT
         ENDIF
      ENDIF               
      write(nuout,*) 'Reading climatological SST for time ',sstclim_time
      
      start(3)=NINT((sstclim_time-first_timein)*kpp_const_fields%spd/
     +     (kpp_const_fields%dto*ndtupdsst))+1
c      write(nuout,*) 'SSTs are being read from position',start(3)
      status=NF_GET_VAR1_REAL(ncid,time_varid,start(3),time_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      IF (abs(time_in-sstclim_time) .GT. 
     +     0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
         write(nuerr,*) 'Cannot find time,',sstclim_time,
     &        'in SST climatology file'
         write(nuerr,*) 'The closest I came was',time_in
         CALL MIXED_ABORT
      ENDIF      

      status=NF_GET_VARA_REAL(ncid,varid,start,count
     $     ,var_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
c      WRITE(nuout,*) 'SSTs have been read from position',start(3)
      status=NF_CLOSE(ncid)
c
c     KPP expects temperatures in CELSIUS.  If climatological SSTs are
c     in Kelvin, subtract 273.15.
c
      offset_sst=0.
      DO ix=1,sst_nx
         DO iy=1,sst_ny
            IF (var_in(ix,iy,1) .gt. 200 .and. 
     &           var_in(ix,iy,1) .lt. 400)
     &           offset_sst = 273.15
         END DO
      END DO
      
      DO ix=1,sst_nx
         DO iy=1,sst_ny

c     This is for the aqua-planet version; it sets the first
c     two rows as land points.  Comment out if not doing aqua-planet.
c     IF ((iy .LE. 2) .OR. (iy .GE. NY_GLOBE-1)) THEN
c     sst_in(ix,iy,1)=-1073741824.0
c     ice_in(ix,iy)=sst_in(ix,iy,1)
c     usf_in(ix,iy)=0.0
c     vsf_in(ix,iy)=0.0
c     ELSE
c     sst_in(ix,iy,1)=var_in(ix,iy,1)            
            sst_in(ix,iy,1)=var_in(ix,iy,1)-offset_sst          
            IF (.NOT. L_CLIMICE) ice_in(ix,iy,1)=0.0
            IF (.NOT. L_CLIMCURR) THEN
               usf_in(ix,iy)=0.0
               vsf_in(ix,iy)=0.0
            ENDIF
c     ENDIF
         ENDDO
      ENDDO

      WRITE(nuout,*) 'Finished read_sstin'

      END

      SUBROUTINE read_icein(kpp_3d_fields,kpp_const_fields)
      
c     Read in ice concentrations from a user-provided netCDF file.
c     Called _only_ if L_CLIMICE is TRUE in 3D_ocn.nml.
c     Called every time the climatological SST is updated.
c     Probably only necessary when coupling the model to an atmosphere
c     that requires realistic ice concentrations.
c     Written by Nick Klingaman, 11/01/08.
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
      include 'constants.com'
      include 'couple.com'
      include 'times.com'
      include 'timocn.com'
      include 'sstclim.com'      
c      include 'location.com'
      
      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      INTEGER ice_nx,ice_ny
#ifdef COUPLE
      PARAMETER(ice_nx=NX_GLOBE,ice_ny=NY_GLOBE)
#else
      PARAMETER(ice_nx=NX,ice_ny=NY)
#endif     
      REAL sst_in(ice_nx,ice_ny,1),ice_in(ice_nx,ice_ny,1),
     &     icedepth_in(ice_nx,ice_ny,1),snowdepth_in(ice_nx,ice_ny,1),
     &     usf_in(ice_nx,ice_ny),vsf_in(ice_nx,ice_ny),
     &     max_ice,min_ice
      REAL*4 var_in(ice_nx,ice_ny,1),iceclim_time,first_timein,
     &     last_timein,time_in,latitudes(NY_GLOBE),longitudes(NX_GLOBE)
      INTEGER count(3),start(3)
      INTEGER ix,iy,status,ncid,varid,time_varid,lon_varid,lat_varid,
     &     time_dimid,lon_dimid,lat_dimid,ntime_file,nlon_file,
     &     nlat_file
      CHARACTER(LEN=30) tmp_name

      COMMON /save_sstin/ sst_in,ice_in,icedepth_in,snowdepth_in,
     &     usf_in,vsf_in
      
c     Set start and count to read a global field if coupled,
c     or a regional field if not coupled.
      count(1)=ice_nx
      count(2)=ice_ny
      count(3)=1
      start(1)=1
      start(2)=1
      start(3)=1

c     Open the netCDF file and find the correct latitude,
c     longitude and time.
      WRITE(nuout,*) 'Opening ice input file'
      status=NF_OPEN(icein_file,0,ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      WRITE(nuout,*) 'In read_icein'

#ifndef COUPLE
      CALL determine_netcdf_boundaries(ncid,'ice climatology',
     &     'latitude','longitude','time',kpp_3d_fields%dlon(1),
     +     kpp_3d_fields%dlat(1),start(1),
     &     start(2),first_timein,last_timein,time_varid)
#else
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid,'time',time_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIMID(ncid,'time',time_dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIM(ncid,time_dimid,tmp_name,ntime_file)
      status=NF_GET_VAR1_REAL(ncid,
     &     time_varid,start(3),first_timein)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR1_REAL(ncid,
     &     time_varid,ntime_file,last_timein)
#endif      

      status=NF_INQ_VARID(ncid,'iceconc',varid)
      iceclim_time=kpp_const_fields%time+0.5*kpp_const_fields%dto/
     +     kpp_const_fields%spd*ndtupdice     
     
      IF (iceclim_time .gt. last_timein) THEN
         IF (L_PERIODIC_CLIMICE) THEN 
            DO WHILE (iceclim_time .gt. last_timein)
               iceclim_time=iceclim_time-climice_period
            ENDDO
         ELSE
            WRITE(nuout,*) 'Time for which to read ice exceeds
     &the last time in the netCDF file and L_PERIODIC_CLIMICE has
     &not been specified.  Attempting to read ice will lead to
     &an error, so aborting now ...'
            CALL MIXED_ABORT
         ENDIF
      ENDIF               
      
      write(nuout,*) 'Reading climatological ICECONC for time ',
     &     iceclim_time
      start(3)=NINT((iceclim_time-first_timein)*kpp_const_fields%spd/
     +     (kpp_const_fields%dto*ndtupdice))+1      
      status=NF_GET_VAR1_REAL(ncid,time_varid,start(3),time_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (abs(time_in-iceclim_time) .GT. 
     +     0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
         write(nuerr,*) 'Cannot find time,',iceclim_time,
     &        'in ice concentration climatology file'
         write(nuerr,*) 'The closest I came was',time_in
         CALL MIXED_ABORT
      ENDIF
      write(nuout,*) 'Ice concentrations are being read from position',
     &     start(3)
      WRITE(nuout,*) 'Start = ',start,'Count = ',count
      status=NF_GET_VARA_REAL(ncid,varid,start,count
     &     ,var_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Ice concentrations have been read from position',
     &     start(3)

      max_ice = -1000.
      min_ice = 1000.
      DO ix=1,ice_nx
         DO iy=1,ice_ny
            ice_in(ix,iy,1) = var_in(ix,iy,1)
            IF (ice_in(ix,iy,1) .gt. max_ice) max_ice = ice_in(ix,iy,1)
            IF (ice_in(ix,iy,1) .lt. min_ice) min_ice = ice_in(ix,iy,1)
         ENDDO
      ENDDO

      IF (L_CLIM_ICE_DEPTH) THEN 
         status=NF_INQ_VARID(ncid,'icedepth',varid)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         WRITE(nuout,*) 'Reading climatological ICEDEPTH for time ',
     &        iceclim_time         
         WRITE(nuout,*) 'Ice depths are being read from ',
     &        'position', start(3)
         WRITE(nuout,*) 'Start = ',start,'Count = ',count
         status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)      

         DO ix=1,ice_nx
            DO iy=1,ice_ny
               icedepth_in(ix,iy,1)=var_in(ix,iy,1)
            ENDDO
         ENDDO
      ENDIF

      IF (L_CLIM_SNOW_ON_ICE) THEN
         status=NF_INQ_VARID(ncid,'snowdepth',varid)
         IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
         WRITE(nuout,*) 'Reading climatological SNOWDEPTH for time ',
     &        iceclim_time
         WRITE(nuout,*) 'Snow depths on sea ice are being read from ',
     &        'position', start(3)
         WRITE(nuout,*) 'Start = ',start,'Count = ',count
         status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
         IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
         
         DO ix=1,ice_nx
            DO iy=1,ice_ny
               snowdepth_in(ix,iy,1)=var_in(ix,iy,1)
            ENDDO
         ENDDO
      ENDIF
      
      status=NF_CLOSE(ncid)

      RETURN
      END

      SUBROUTINE read_surface_currents(kpp_3d_fields,kpp_const_fields)

      IMPLICIT NONE

      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <netcdf.inc>
! Automatically includes paramter.inc!
#include <kpp_3d_type.com>
#include <constants.com>
#include <couple.com>
#include <times.com>
#include <timocn.com>
c#include <location.com>
#include <sstclim.com>
#include <currclim.com>

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      
      INTEGER curr_nx,curr_ny
#ifdef COUPLE
      PARAMETER(curr_nx=NX_GLOBE,curr_ny=NY_GLOBE)
#else
      PARAMETER(curr_nx=NX,curr_ny=NY)
#endif     
      REAL sst_in(curr_nx,curr_ny,1),ice_in(curr_nx,curr_ny,1),
     &     icedepth_in(curr_nx,curr_ny,1),
     &     snowdepth_in(curr_nx,curr_ny,1),
     &     usf_in(curr_nx,curr_ny),vsf_in(curr_nx,curr_ny),
     &     max_ice,min_ice
      REAL*4 var_in(curr_nx,curr_ny),currclim_time,first_timein,
     &     time_in,latitudes(NY_GLOBE),longitudes(NX_GLOBE),last_timein
      INTEGER count(3),start(3)
      INTEGER ix,iy,status,ncid,varid,time_varid,time_dimid,ntime_file
      CHARACTER(LEN=30) tmp_name

      COMMON /save_sstin/ sst_in,ice_in,icedepth_in,snowdepth_in,
     &     usf_in,vsf_in
            
c     Set start and count to read a global field if coupled,
c     or a regional field if not coupled.
      count(1)=curr_nx
      count(2)=curr_ny
      count(3)=1
      start(1)=1
      start(2)=1
      start(3)=1

c     Open the netCDF file and find the correct time.
      status=NF_OPEN(currin_file,0,ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)


#ifndef COUPLE
      CALL determine_netcdf_boundaries(ncid,'current climatology',
     &     'latitude','longitude','time',dlon(1),dlat(1),start(1),
     &     start(2),first_timein,last_timein,time_varid)
#else
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_VARID(ncid,'time',time_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIMID(ncid,'time',time_dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIM(ncid,time_dimid,tmp_name,ntime_file)
      status=NF_GET_VAR1_REAL(ncid,
     &     time_varid,start(3),first_timein)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR1_REAL(ncid,
     &     time_varid,ntime_file,last_timein)
#endif      

      currclim_time=time+0.5*dto/spd*ndtupdsst
      write(nuout,*) 'Reading climatological USF for time ',
     &     currclim_time
      start(3)=NINT((currclim_time-first_timein)*spd/(dto*ndtupdcurr))+1      
      status=NF_GET_VAR1_REAL(ncid,time_varid,start(3),time_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (abs(time_in-currclim_time) .GT. 0.01*kpp_const_fields%dtsec/
     +     kpp_const_fields%spd) THEN
         write(nuerr,*) 'Cannot find time,',currclim_time,
     &        'in curr concentration climatology file'
         write(nuerr,*) 'The closest I came was',time_in
         CALL MIXED_ABORT
      ENDIF
      write(nuout,*) 'Zonal currents are being read from position',
     &     start(3)
      WRITE(nuout,*) 'Start = ',start,'Count = ',count
      status=NF_GET_VARA_REAL(ncid,varid,start,count
     &     ,var_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Zonal currents have been read from position',
     &     start(3)

      DO ix=1,curr_nx
         DO iy=1,curr_ny
            usf_in(ix,iy)=var_in(ix,iy)
         ENDDO
      ENDDO

      status=NF_INQ_VARID(ncid,'vcurr',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Reading climatological VSF for time ',
     &     currclim_time
      WRITE(nuout,*) 'Meridional currents are being read from position',
     &     start(3)
      status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Meridional currents have been read from position',
     &     start(3)

      DO ix=1,curr_nx
         DO iy=1,curr_ny
            vsf_in(ix,iy)=var_in(ix,iy)
         ENDDO
      ENDDO
  
      status=NF_CLOSE(ncid)

      RETURN
      END

      SUBROUTINE read_bottom_temp(kpp_3d_fields,kpp_const_fields,
     +     bottom_temp)

      IMPLICIT NONE

      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
#include <constants.com>
#include <couple.com>
#include <times.com>
#include <timocn.com>
#include <bottomclim.com>
c#include <location.com>

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      INTEGER status,ncid
      REAL bottom_temp(NPTS),offset_temp
      REAL*4 var_in(NX,NY,1),time_in,
     & first_timein, bottomclim_time,latitudes(NY_GLOBE), 
     & longitudes(NX_GLOBE), last_timein
      INTEGER varid,time_varid,lat_varid,lon_varid,time_dimid,lat_dimid,
     &     lon_dimid,nlat_file,nlon_file,ntime_file

      INTEGER count(3),start(3)
      INTEGER ix,iy,ipoint
      CHARACTER(LEN=30) tmp_name

c     Set start and count to read a regional field.
      count(1)=NX
      count(2)=NY
      count(3)=1
      start(1)=1
      start(2)=1
      start(3)=1

c     Open the netCDF file and find the latitude and longitude 
c     boundaries in the input file.
      status=NF_OPEN(bottomin_file,0,ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      status=NF_INQ_VARID(ncid,'T',varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      CALL determine_netcdf_boundaries(ncid,'bottom temp climatology',
     &     'latitude','longitude','time',kpp_3d_fields%dlon(1),
     +     kpp_3d_fields%dlat(1),start(1),
     &     start(2),first_timein,last_timein,time_varid)

      bottomclim_time=kpp_const_fields%time+0.5*kpp_const_fields%dto/
     +     kpp_const_fields%spd*ndtupdbottom
      IF (bottomclim_time .gt. last_timein) THEN
         IF (L_PERIODIC_BOTTOM_TEMP) THEN 
            DO WHILE (bottomclim_time .gt. last_timein)
               bottomclim_time=bottomclim_time-bottom_temp_period
            ENDDO
         ELSE
            WRITE(nuout,*) 'Time for which to read bottom temperature
     &exceeds the last time in the netCDF file and
     &L_PERIODIC_BOTTOM_TEMP has not been specified.  Attempting to
     &read bottom temperature will lead to an error, so aborting
     &now...'
            CALL MIXED_ABORT
         ENDIF
      ENDIF               
 
      write(nuout,*) 'Reading climatological bottom temp for time ',
     &     bottomclim_time
      start(3)=NINT((bottomclim_time-first_timein)*kpp_const_fields%spd/
     &     (kpp_const_fields%dto*ndtupdbottom))+1
      write(nuout,*) 'Bottom temperatures are being read from position',
     &     start(3)

      status=NF_GET_VAR1_REAL(ncid,time_varid,start(3),time_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (abs(time_in-bottomclim_time) .GT. 0.01*kpp_const_fields%dtsec/
     +     spd) THEN
         write(nuerr,*) 'Cannot find time',bottomclim_time,
     &        'in bottom temperature climatology file'
         write(nuerr,*) 'The closest I came was',time_in
         CALL MIXED_ABORT
      ENDIF
      status=NF_GET_VARA_REAL(ncid,varid,start,count
     &     ,var_in)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Bottom temperatures have been read from position',
     &     start(3)
      WRITE(nuout,*) 'First column of bottom temperatures: ',
     +     var_in(1,:,1)
c
c     KPP expects temperatures in CELSIUS.  If climatological bottom 
c     temperatures are in Kelvin, subtract 273.15.
c     
      offset_temp=0.
      ix=1
      iy=1
      DO WHILE (offset_temp.EQ.0.AND.ix.LE.NX)
         DO iy=1,NY
            IF (var_in(ix,iy,1) .gt. 200 .and. 
     &           var_in(ix,iy,1) .lt. 400)
     &           offset_temp = 273.15         
         END DO
         ix=ix+1
      ENDDO

      DO ix=1,NX
         DO iy=1,NY
            ipoint=(iy-1)*NX+ix
            bottom_temp(ipoint) = var_in(ix,iy,1)
            bottom_temp(ipoint) = bottom_temp(ipoint)-offset_temp
         ENDDO
      ENDDO     

      status=NF_CLOSE(ncid)

      RETURN
      END

      SUBROUTINE read_salinity(kpp_3d_fields,kpp_const_fields)
      
      IMPLICIT NONE
      INTEGER nuout,nuerr,start(4),count(4)
      INTEGER ix,iy,iz,ipoint,sal_varid,status,lat_varid,lon_varid,
     +     z_varid,z_dimid,time_varid,sal_ncid,k,lat_dimid,lon_dimid,
     +     time_dimid,nlon_file,nlat_file,ntime_file,nz_file

      PARAMETER (nuout=6,nuerr=0)
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
c#include <location.com>
#include <relax_3d.com>
#include <times.com>
#include <timocn.com>
#include <constants.com>
#include <ocn_paras.com>
      
      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      REAL*4 ixx,jyy,first_timein,time_in,
     +     sal_time,ndays_upd_sal,last_timein
      CHARACTER(LEN=30) tmp_name
      REAL*4, allocatable :: sal_in(:,:,:,:),latitudes(:),longitudes(:),
     +     z(:)
c      COMMON /save_sal/ sal_clim
c
c     Read in a NetCDF file containing a 
c     time-varying salinity field at every model vertical level.
c     Frequency of read is controlled by ndtupdsal in the namelist
c
c     NPK 12/02/08
c      
      allocate(sal_in(NX,NY,NZP1,1))
      allocate(longitudes(NX_GLOBE))
      allocate(latitudes(NY_GLOBE))
      allocate(z(NZP1))

      WRITE(nuout,*) 'Trying to open salinity input file ',sal_file
      status=NF_OPEN(sal_file,0,sal_ncid)
      IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Opened salinity input file ',sal_ncid
      
      count=(/NX,NY,NZP1,1/)
      start=(/1,1,1,1/)

      status=NF_INQ_VARID(sal_ncid,'z',z_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)  
      status=NF_INQ_DIMID(sal_ncid,'z',z_dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIM(sal_ncid,z_dimid,tmp_name,nz_file)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (NZP1.ne.nz_file) THEN
         WRITE(nuout,*) 'Input file for salinity climatology does ',
     &        'not have the correct number of vertical levels. ',
     &        'It should have ',NZP1,' but instead has ',nz_file
         CALL MIXED_ABORT
      ELSE
         status=NF_GET_VAR_REAL(sal_ncid,z_varid,z)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         WRITE(nuout,*) 'Read in depths from the salinity climatology ',
     &        'input file'
      ENDIF

      WRITE(6,*) kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1)
      CALL determine_netcdf_boundaries(sal_ncid,'salinity clim',
     &     'latitude','longitude','t',kpp_3d_fields%dlon(1),
     +     kpp_3d_fields%dlat(1),start(1),
     &     start(2),first_timein,last_timein,time_varid)

      status=NF_INQ_VARID(sal_ncid,'salinity',sal_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      ndays_upd_sal = ndtupdsal*kpp_const_fields%dto/
     +     kpp_const_fields%spd
c      WRITE(nuout,*) ndays_upd_sal,FLOOR(time)*NINT(kpp_const_fields%spd),
c     &     ndtupdsal*NINT(dto),0.5*dto/spd*ndtupdsal
      sal_time=(ndays_upd_sal)*
     &     (FLOOR(kpp_const_fields%time)*NINT(kpp_const_fields%spd)/
     +     (ndtupdsal*NINT(kpp_const_fields%dto)))+
     &     (0.5*kpp_const_fields%dto/kpp_const_fields%spd*ndtupdsal)
      WRITE(nuout,*) sal_time,last_timein
      
      IF (sal_time .gt. last_timein) THEN
         IF (L_PERIODIC_SAL) THEN 
            DO WHILE (sal_time .gt. last_timein)
               sal_time=sal_time-sal_period
            ENDDO
         ELSE
            WRITE(nuout,*) 'Time for which to read the salinity
     & climatology exceeds the last time in the netCDF file
     & and L_PERIODIC_SAL has not been specified.
     & Attempting to read salinity climatology will lead to an error, so
     & aborting now ...'
            CALL MIXED_ABORT
         ENDIF
      ENDIF               

      write(nuout,*) 'Reading salinity for time ',sal_time
      start(4)=NINT((sal_time-first_timein)*kpp_const_fields%spd/
     +     (kpp_const_fields%dto*ndtupdsal))+1
      write(nuout,*) 'Salinity values are being read from position',
     &     start(4)
      status=NF_GET_VAR1_REAL(sal_ncid,time_varid,start(4),time_in)
      
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (abs(time_in-sal_time) .GT. 0.01*kpp_const_fields%dtsec/
     +     kpp_const_fields%spd) THEN
         write(nuerr,*) 'Cannot find time',sal_time,
     &        'in salinity climatology input file'
         write(nuerr,*) 'The closest I came was',time_in
         CALL MIXED_ABORT
      ENDIF
      status=NF_GET_VARA_REAL(sal_ncid,sal_varid,start,count
     &     ,sal_in)
      write(nuout,*) 'Salinity climatology data have been read from '//
     &     'position',start(4)
      
c
c     Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
c     into one long array with dimension NPTS.
c         
      DO ix=1,NX
         DO iy=1,NY
            ipoint=(iy-1)*nx+ix
            DO k=1,NZP1
c
c     Subtract reference salinity from climatology, for compatability with
c     salinity values stored in main model.
c               
               kpp_3d_fields%sal_clim(ipoint,k)=sal_in(ix,iy,k,1)-
     +              kpp_3d_fields%Sref(ipoint)
            ENDDO
         ENDDO
      ENDDO

      status=NF_CLOSE(sal_ncid)
      deallocate(sal_in)
      deallocate(longitudes)
      deallocate(latitudes)
      deallocate(z)

      RETURN
      END

      SUBROUTINE read_ocean_temperatures(kpp_3d_fields,kpp_const_fields)
      IMPLICIT NONE

#include <kpp_3d_type.com>
#include <netcdf.inc>
#include <relax_3d.com>
      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      INTEGER nuout,nuerr,start(4),count(4)
      INTEGER ix,iy,iz,ipoint,ocnT_varid,status,lat_varid,lon_varid,
     +     z_varid,z_dimid,time_varid,ocnT_ncid,k,lat_dimid,lon_dimid,
     +     time_dimid,nlon_file,nlat_file,ntime_file,nz_file

      PARAMETER (nuout=6,nuerr=0)
      REAL*4, allocatable :: ocnT_in(:,:,:,:),latitudes(:),
     +     longitudes(:),z(:)
      REAL*4 ixx,jyy,first_timein,time_in,
     +     ocnT_time,ndays_upd_ocnT,last_timein
      CHARACTER(LEN=30) tmp_name

      allocate(ocnT_in(NX,NY,NZP1,1))
      allocate(longitudes(NX_GLOBE))
      allocate(latitudes(NY_GLOBE))
      allocate(z(NZP1))

      WRITE(nuout,*) 'Trying to open ocean temperature input file ',
     +     ocnT_file
      status=NF_OPEN(ocnT_file,0,ocnT_ncid)
      IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Opened ocean temperature input file ',ocnT_ncid
      
      count=(/NX,NY,NZP1,1/)
      start=(/1,1,1,1/)

      status=NF_INQ_VARID(ocnT_ncid,'z',z_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)  
      status=NF_INQ_DIMID(ocnT_ncid,'z',z_dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIM(ocnT_ncid,z_dimid,tmp_name,nz_file)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (NZP1.ne.nz_file) THEN
         WRITE(nuout,*) 'Input file for ocean temperature climatology ',
     &        'does not have the correct number of vertical levels. ',
     &        'It should have ',NZP1,' but instead has ',nz_file
         CALL MIXED_ABORT
      ELSE
         status=NF_GET_VAR_REAL(ocnT_ncid,z_varid,z)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         WRITE(nuout,*) 'Read in depths from the ocean temperature ',
     &        'climatology input file'
      ENDIF

      WRITE(6,*) kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1)
      CALL determine_netcdf_boundaries(ocnT_ncid,'ocean temp clim',
     &     'latitude','longitude','t',kpp_3d_fields%dlon(1),
     +     kpp_3d_fields%dlat(1),start(1),
     &     start(2),first_timein,last_timein,time_varid)

      status=NF_INQ_VARID(ocnT_ncid,'temperature',ocnT_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      ndays_upd_ocnT = ndtupdocnT*kpp_const_fields%dto/
     +     kpp_const_fields%spd
c      WRITE(nuout,*) ndays_upd_ocnT,FLOOR(time)*NINT(kpp_const_fields%spd),
c     &     ndtupdocnT*NINT(dto),0.5*dto/spd*ndtupdocnT
      ocnT_time=(ndays_upd_ocnT)*
     &     (FLOOR(kpp_const_fields%time)*NINT(kpp_const_fields%spd)/
     +     (ndtupdocnT*NINT(kpp_const_fields%dto)))+
     &     (0.5*kpp_const_fields%dto/kpp_const_fields%spd*ndtupdocnT)
      WRITE(nuout,*) ocnT_time,last_timein
      
      IF (ocnT_time .gt. last_timein) THEN
         IF (L_PERIODIC_OCNT) THEN 
            DO WHILE (ocnT_time .gt. last_timein)
               ocnT_time=ocnT_time-ocnT_period
            ENDDO
         ELSE
            WRITE(nuout,*) 'Time for which to read the ocean
     & temperatures exceeds the last time in the netCDF file
     & and L_PERIODIC_OCNT has not been specified.
     & Attempting to read ocean temperatures will lead to an error, so
     & aborting now ...'
            CALL MIXED_ABORT
         ENDIF
      ENDIF               

      write(nuout,*) 'Reading ocean temperature for time ',ocnT_time
      start(4)=NINT((ocnT_time-first_timein)*kpp_const_fields%spd/
     +     (kpp_const_fields%dto*ndtupdocnT))+1
      write(nuout,*) 'Ocean temperature values are being read from ',
     &     'position',start(4)
      status=NF_GET_VAR1_REAL(ocnT_ncid,time_varid,start(4),time_in)
      
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      IF (abs(time_in-ocnT_time) .GT. 0.01*kpp_const_fields%dtsec/
     +     kpp_const_fields%spd) THEN
         write(nuerr,*) 'Cannot find time',ocnT_time,
     &        'in ocean temperature climatology input file'
         write(nuerr,*) 'The closest I came was',time_in
         CALL MIXED_ABORT
      ENDIF
      status=NF_GET_VARA_REAL(ocnT_ncid,ocnT_varid,start,count
     &     ,ocnT_in)
      write(nuout,*) 'Ocean temperature climatology data have been '//
     &     'read from position',start(4)
      
c
c     Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
c     into one long array with dimension NPTS.
c         
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
      END

      SUBROUTINE determine_netcdf_boundaries(ncid,file_description,
     &     latitude_name,longitude_name,time_name,start_lon,start_lat,
     &     offset_lon,offset_lat,first_time,last_time,time_varid)

      IMPLICIT NONE
      
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <parameter.inc>
#include <netcdf.inc>

      INTEGER ncid,offset_lon,offset_lat,lon_dimid,lon_varid,
     &     lat_dimid,lat_varid,time_dimid,time_varid,ix,iy
      REAL start_lon,start_lat,first_time,last_time      
      CHARACTER(*) file_description,latitude_name,longitude_name,
     &     time_name
      CHARACTER(LEN=30) tmp_name

      INTEGER nlat_file,nlon_file,ntime_file,status
      REAL*4 longitudes(NX_GLOBE),latitudes(NY_GLOBE)
      
      status=NF_INQ_DIMID(ncid,latitude_name,lat_dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIM(ncid,lat_dimid,tmp_name,nlat_file)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)      
      status=NF_INQ_VARID(ncid,latitude_name,lat_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Inquired on '//latitude_name,lat_varid

      status=NF_INQ_DIMID(ncid,longitude_name,lon_dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIM(ncid,lon_dimid,tmp_name,nlon_file)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)      
      status=NF_INQ_VARID(ncid,longitude_name,lon_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Inquired on '//longitude_name,lon_varid

      status=NF_GET_VAR_REAL(ncid,lat_varid,latitudes(1:nlat_file))
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Read ',nlat_file,' latitudes'
      status=NF_GET_VAR_REAL(ncid,lon_varid,longitudes(1:nlon_file))
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Read ',nlon_file,' longitudes'

      ix=1
      DO WHILE (abs(longitudes(ix)-start_lon) .GT. 1.e-3)
         IF (ix .gt. nlon_file) THEN
            WRITE(nuout,*) 'Could not find starting ',
     &           'longitude in the '//file_description//' file'
            CALL MIXED_ABORT
         ENDIF
         ix=ix+1
      ENDDO
      offset_lon=ix
       
      iy=1
      DO WHILE (abs(latitudes(iy)-start_lat) .GT. 1.e-3)
         IF (iy .gt. nlat_file) THEN
            WRITE(nuout,*) 'Could not find starting ',
     &           'latitude in the '//file_description//' file'
            CALL MIXED_ABORT
         ENDIF
         iy=iy+1
      ENDDO
      offset_lat=iy

c     Find the first time and last time
      status=NF_INQ_VARID(ncid,time_name,time_varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR1_REAL(ncid,
     &     time_varid,1,first_time)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIMID(ncid,time_name,time_dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_INQ_DIM(ncid,time_dimid,tmp_name,ntime_file)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_GET_VAR1_REAL(ncid,
     &     time_varid,ntime_file,last_time)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      RETURN
      END
