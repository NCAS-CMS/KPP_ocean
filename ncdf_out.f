      SUBROUTINE init_output(filename,ncid,kpp_3d_fields,
     +     kpp_const_fields)
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
c#include <landsea.com>
      include 'vert_pgrid.com'
      include 'output.com'
c      include 'location.com'

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      INTEGER k,ix,iy,ncid
      CHARACTER*11 varname(N_VAROUTS),singname(N_SINGOUTS)
      CHARACTER*50 longname(N_VAROUTS),singlong(N_SINGOUTS),
     &     filename
      CHARACTER*15 units(N_VAROUTS) ,singunits(N_SINGOUTS)
      CHARACTER*6 type

      integer vert_dimid(N_VAROUTS)
      integer dims(4)
      integer zflag(N_VAROUTS) ! 1 for z-levels, 0 for d-levels 

      INTEGER status

      REAL*4 ZOUT(NZP1),ALON(NX),ALAT(NY)
      REAL*4 delta
      
      DATA varname/'u','v','T','S','B',
     &     'wu','wv','wT','wS','wB',
     &     'wTnt',
     &     'difm','dift','difs',
     &     'rho','cp','scorr','Rig','dbloc','Shsq','tinc_fcorr',
     &     'fcorr_z','sinc_fcorr'/
      DATA zflag/1,1,1,1,1,
     &     0,0,0,0,0,
     &     0,
     &     0,0,0,
     &     1,1,1,1,1,1,1,
     &     1,1/
      DATA longname/
     &     'Zonal velocity',
     &     'Meridional velocity',
     &     'Temperature',
     &     'Salinity',
     &     'Buoyancy',
     &     'Turbulent Zonal Velocity Flux',
     &     'Turbulent Meridional Velocity Flux',
     &     'Turbulent Temperature Flux',
     &     'Turbulent Salinity Flux',
     &     'Turbulent Buoyancy Flux',
     &     'Non-Turbulent Temperature Flux',
     &     'Diffusion Coefficient (Momentum)',
     &     'Diffusion Coefficient (Temperature)',
     &     'Diffusion Coefficient (Salinity)',
     &     'Density',
     &     'Specific Heat Capacity',
     &     'Salinity correction (with depth)',
     &     'Local Richardson Number in kpp.f',
     &     'Local delta buoyancy in kpp.f',
     &     'Local shear-squared term in kpp.f',
     &     'Temperature increment flux correction',
     &     'Heat correction as flux (dT/dt*rho*cp)',
     &     'Salinity increment flux correction'/
      DATA units/
     &     'm/s',
     &     'm/s',
     &     'degC',
     &     'o/oo',
     &     'm/s^2',
     &     'm^2/s^2',
     &     'm^2/s^2',
     &     'degC m/s',
     &     'o/oo m/s',
     &     'm^2/s^3',
     &     'degC m/s',
     &     'm^2/s',
     &     'm^2/s',
     &     'm^2/s',
     &     'kg/m^3',
     &     'J/kg/K',
     &     'o/oo/s',
     &     'unitless',
     &     'm/s^2',
     &     'm^2/s^2',
     &     'K/timestep',
     &     'W/m^3',
     &     'o/oo/timestep'/
      DATA singname/
     &     'hmix',
     &     'fcorr',
     &     'taux_in',
     &     'tauy_in',
     &     'solar_in',
     &     'nsolar_in',
     &     'PminusE_in',
     &     'cplwght', 
     &     'freeze_flag',
     &     'comp_flag'/
#ifdef COUPLE
#ifdef OASIS2
      DATA type /'OASIS2'/
#else
#ifdef OASIS3
      DATA type /'OASIS3'/
#else
#ifdef GFS
      DATA type /' GFS  '/
#endif /*GFS*/
#endif /*OASIS3*/
#endif /*OASIS2*/
#else
      DATA type /'netCDF'/
#endif /*COUPLE*/

      DATA singlong/
     &     'Mixed Layer Depth',
     &     'Flux Correction',
     &     'Zonal wind stress from ',
     &     'Meridional wind stress from ',
     &     'Solar from ',
     &     'Non-solar from ',
     &     'P minus E from ',
     &     'Coupling weight',
     &     'Fraction of levels below freezing',
     &     'Number of integrations (<0 = isothermal reset)'/

      DATA singunits/
     &     'm',
     &     'W/m^2',
     &     'N/m^2',
     &     'N/m^2',
     &     'W/m^2',
     &     'W/m^2',
     &     'mm/s',
     &     'none',
     &     'fraction',
     &     'unitless'/
      
      nout=1
      singlong(3:7)=singlong(3:7)//type

      DO ix=1,nx
         alon(ix)=kpp_3d_fields%dlon(ix)
      ENDDO
      DO iy=1,ny
         alat(iy)=kpp_3d_fields%dlat((iy-1)*nx+1)
      ENDDO

      status=NF_CREATE(filename, nf_clobber, ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Output file ',filename,' created successfully.'

      delta=0.0
      IF (NX .GT. 1 ) delta=alon(2)-alon(1)
      CALL MY_NCDF_DEF_DIM (
     &     ncid,londim,nx,lon_id,'longitude','deg',delta,' ')
c      WRITE(nuout,*) 'X dimension created successfully for ',filename
      delta=0.0
      IF (NY .GT. 1) delta=alat(2)-alat(1)
      CALL MY_NCDF_DEF_DIM (
     &     ncid,latdim,ny,lat_id,'latitude','deg',delta,' ')
c      WRITE(nuout,*) 'Y dimension created successfully for ',filename
      CALL MY_NCDF_DEF_DIM (
     &     ncid,zdim,NZP1,z_id,'z','m',0.0,' ')
c      WRITE(nuout,*) 'Z dimension created successfully for ',filename
      CALL MY_NCDF_DEF_DIM (
     &     ncid,hdim,NZP1,h_id,'h','m',0.0,'Layer Thickness')
c      WRITE(nuout,*) 'H dimension created successfully for ',filename
      CALL MY_NCDF_DEF_DIM (
     &     ncid,ddim,NZP1,d_id,'d','m',0.0,'Depth of Interfaces')
c      WRITE(nuout,*) 'D dimension created successfully for ',filename
      delta=dtout
      CALL MY_NCDF_DEF_DIM (
     &     ncid,timdim,NF_UNLIMITED,time_id,'time',
     &     'days',delta,' ')
c      WRITE(nuout,*) 'T dimension created successfully for ',filename
      
      DO k=1,N_VAROUTS
         IF (zflag(k) .EQ. 1) THEN 
            vert_dimid(k)=zdim
         ELSEIF (zflag(k) .EQ. 0) THEN
            vert_dimid(k)=ddim
         ELSE
            write(nuerr,*) 'Incorrect value for zflag, for ',
     &           longname(k)
            CALL MIXED_ABORT
         ENDIF
      ENDDO
      
      dims(1)=londim
      dims(2)=latdim
      dims(4)=timdim
      DO k=1,N_VAROUTS
         dims(3)=vert_dimid(k)
         IF (L_VAROUT(k)) THEN
            CALL MY_NCDF_DEF_VAR (
     &           ncid,varid(k),4,dims,varname(k),units(k),
     &           longname(k))
c            WRITE(nuout,*) 'Vector',k,'created successfully for ',
c     &           filename,'with varid=',varid(k)
         ENDIF
      ENDDO
      
      dims(3)=timdim
      DO k=1,N_SINGOUTS
         IF (L_SINGOUT(k)) THEN
            CALL MY_NCDF_DEF_VAR (
     &           ncid,singid(k),3,dims,singname(k),singunits(k),
     &           singlong(k))
c     WRITE(nuout,*) 'Scalar ',k,' created successfully for ',
c     &           filename,'with varid=',singid(k)
         ENDIF
      ENDDO
      
      status=NF_ENDDEF(ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      status=NF_PUT_VAR_REAL(ncid,lon_id,alon)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      status=NF_PUT_VAR_REAL(ncid,lat_id,alat)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      DO k=1,NZP1
         ZOUT(k)=kpp_const_fields%zm(k)
      ENDDO
      status=NF_PUT_VAR_REAL(ncid,z_id,ZOUT)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      DO k=1,NZP1
         ZOUT(k)=kpp_const_fields%hm(k)
      ENDDO
      status=NF_PUT_VAR_REAL(ncid,h_id,ZOUT)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      DO k=1,NZP1
         ZOUT(k)=kpp_const_fields%dm(k-1)
      ENDDO
      status=NF_PUT_VAR_REAL(ncid,d_id,ZOUT)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      write(nuout,*) 'NCDF output file initialized with name ',
     +      output_file
      
c      call output_close
      status=NF_CLOSE(ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      RETURN
      END

      SUBROUTINE output_inst(kpp_3d_fields,kpp_const_fields)
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
      include 'output.com'

      include 'times.com'
      include 'ocn_advec.com'
#include <landsea.com>
#include <relax_3d.com>
      include 'couple.com'
#include <fcorr_in.com>
#include <sfcorr_in.com>

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

      REAL*4, allocatable :: varout(:,:,:), singout(:,:)
      REAL*4 TOUT
      
      INTEGER start(4),count(4),status
      INTEGER k,ivar
      INTEGER ix,iy,ipt
 
      allocate(varout(NX,NY,NZP1))
      allocate(singout(NX,NY))

      count(1)=NX
      count(2)=NY
      count(3)=NZP1
      count(4)=1
      
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=nout

c     NPK 25/2/08.
c     Stop opening and closing output files with each flush of output.
c     Just do it once at the beginning/end of the simulation.
c      call output_open
      TOUT=kpp_const_fields%time
      status=NF_PUT_VAR1_REAL(ncid_out,time_id,nout,TOUT)

c      write(nuout,*) 'Writing output at timestep ',ntime+nstart,
c     +        ' Time=',TOUT

      DO ivar=1,N_VAROUTS
         IF (L_VAROUT(ivar)) THEN
            VAROUT(:,:,:)=missval
            IF (ivar .EQ. 1) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%U(ipt,k,1)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 2) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%U(ipt,k,2)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 3) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%X(ipt,k,1)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 4) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%X(ipt,k,2)+
     +                          kpp_3d_fields%Sref(ipt)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 5) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%buoy(ipt,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 6) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF(kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%wU(ipt,k-1,1)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 7) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF(kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%wU(ipt,k-1,2)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 8) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%wX(ipt,k-1,1)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 9) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%wX(ipt,k-1,2)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 10) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=
     +                          kpp_3d_fields%wX(ipt,k-1,NSP1)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 11) THEN               
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%wXNT(ipt,k-1,1)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 12) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        VAROUT(ix,iy,1)=0.0
                        DO k=2,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%difm(ipt,k-1)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 13) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        VAROUT(ix,iy,1)=0.0
                        DO k=2,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%dift(ipt,k-1)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 14) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        VAROUT(ix,iy,1)=0.0
                        DO k=2,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%difs(ipt,k-1)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 15) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%rho(ipt,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 16) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%cp(ipt,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 17) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%scorr(ipt,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 18) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN 
                        DO k=1,NZ
                           VAROUT(ix,iy,k)=kpp_3d_fields%Rig(ipt,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 19) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZ
                           VAROUT(ix,iy,k)=kpp_3d_fields%dbloc(ipt,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 20) THEN 
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZ
                           VAROUT(ix,iy,k)=kpp_3d_fields%shsq(ipt,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 21) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=
     +                          kpp_3d_fields%Tinc_fcorr(ipt,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO       
            ELSEIF (ivar .EQ. 22) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=kpp_3d_fields%ocnTcorr(ipt,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 23) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=
     +                          kpp_3d_fields%Sinc_fcorr(ipt,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSE
               write(nuerr,*) 'You need to extend the outputs in'//
     &              'output_inst'
            ENDIF
c            WRITE(nuout,*) 'In output_inst for varout, ivar=',ivar
            status=NF_PUT_VARA_REAL(
     &           ncid_out,varid(ivar),start,count,VAROUT)
            IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         ENDIF
      ENDDO
      
      start(3)=start(4)
      count(3)=1
      DO ivar=1,N_SINGOUTS
         IF (L_SINGOUT(ivar)) THEN
            SINGOUT(:,:)=missval
            IF (ivar .EQ. 1) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) 
     &                    SINGOUT(ix,iy)=kpp_3d_fields%hmix(ipt)
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 2) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt))
     &                    SINGOUT(ix,iy)=kpp_3d_fields%fcorr(ipt)                  
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 3) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt))
     &                    SINGOUT(ix,iy)=kpp_3d_fields%sflux(ipt,1,5,0)
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 4) THEN 
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt))
     &                    SINGOUT(ix,iy)=kpp_3d_fields%sflux(ipt,2,5,0)
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 5) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt))
     &                    SINGOUT(ix,iy)=kpp_3d_fields%sflux(ipt,3,5,0)
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 6) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt))
     &                    SINGOUT(ix,iy)=kpp_3d_fields%sflux(ipt,4,5,0)
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 7) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt))
     &                    SINGOUT(ix,iy)=kpp_3d_fields%sflux(ipt,6,5,0)
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 8) THEN
               DO ix=ifirst,ilast
                  DO iy=jfirst,jlast
!                     WRITE(nuout,*) ix,iy,ifirst,ilast,jfirst,jlast
                     ipt=(iy-1)*NX_GLOBE+ix
                     SINGOUT(ix-ifirst+1,iy-jfirst+1)=
     +                    kpp_3d_fields%cplwght(ipt)
                  ENDDO 
               ENDDO
            ELSEIF (ivar .EQ. 9) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*NX+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt))
     &                    SINGOUT(ix,iy)=kpp_3d_fields%freeze_flag(ipt)
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 10) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*NX+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt))
     &                    SINGOUT(ix,iy)=kpp_3d_fields%reset_flag(ipt)
                  ENDDO
               ENDDO
            ELSE
               write(nuerr,*) 'You need to extend the outputs in'//
     &              'output_inst'
            ENDIF
c            WRITE(nuout,*) 'In output_inst for singout, ivar=',ivar
            status=NF_PUT_VARA_REAL(
     &           ncid_out,singid(ivar),start,count,SINGOUT)
            IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         ENDIF
      ENDDO
      start(3)=1

      nout=nout+1
c     NPK 25/2/08
c     Stop opening and closing output files with each flush of output.
c     Just do it once at the beginning/end of the simulation.
c      call output_close
      WRITE(nuout,*) ' Output successfully written'
      RETURN
      END

      SUBROUTINE write_means(kpp_3d_fields,kpp_const_fields,
     +     VEC_mean,SCLR_mean)

      IMPLICIT NONE

      INTEGER nuout,nuerr
      PARAMETER(nuout=6,nuerr=0)

      include 'netcdf.inc'      
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
      include 'output.com'
#include <landsea.com>
      include 'times.com'
      include 'timocn.com'
      include 'ocn_advec.com'
      include 'couple.com'
      include 'constants.com'

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

      REAL :: VEC_mean(NPTS,NZP1,NVEC_MEAN),
     + SCLR_mean(NPTS,NSCLR_MEAN)
      REAL*4,allocatable :: VAROUT(:,:,:), SINGOUT(:,:)
      REAL*4 TOUT

      INTEGER i,ivar,ipt,ix,iy,start(4),count(4),k,status

      allocate(VAROUT(NX,NY,NZP1))
      allocate(SINGOUT(NX,NY))
      TOUT=kpp_const_fields%time-(ndtout_mean*kpp_const_fields%dtsec/
     +     (ndtocn*kpp_const_fields%spd))*0.5
      status=NF_PUT_VAR1_REAL(mean_ncid_out,time_id,nout_mean,TOUT)
c      write(nuout,*) 'Writing mean output at timestep ',ntime+nstart,
c     +     ' Time=',TOUT
      
      count(1)=NX
      count(2)=NY
      count(3)=NZP1
      count(4)=1
      
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=nout_mean
      
      i=1
      DO ivar=1,N_VAROUTS
         IF (L_MEAN_VAROUT(ivar)) THEN
            VAROUT(:,:,:)=missval
            IF (ivar.EQ.1) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar.EQ.2) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 3) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO                     
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 4) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)+
     +                          kpp_3d_fields%Sref(ipt)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 5) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 6) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF(kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
c     Potentially k-1?
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO                    
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 7) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF(kpp_3d_fields%L_OCEAN(ipt)) THEN
c     Potentially k-1?
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO                   
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 8) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
c     Potentially k-1?
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 9) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
c     Potentially k-1?
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 10) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
c     Potentially k-1?
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 11) THEN               
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
c     Potentially k-1?
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 12) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        VAROUT(ix,iy,1)=0.0
                        DO k=2,NZP1
c     Potentially k-1?
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 13) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        VAROUT(ix,iy,1)=0.0
                        DO k=2,NZP1
c     Potentially k-1?
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 14) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN     
                        VAROUT(ix,iy,1)=0.0
                        DO k=2,NZP1
c     Potentially k-1?
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar. EQ. 15) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 16) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO                                       
            ELSEIF (ivar .EQ. 17) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 18) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN 
                        DO k=1,NZ
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 19) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZ
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 20) THEN 
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZ
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 21) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 22) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSEIF (ivar .EQ. 23) THEN
               DO ix=1,nx
                  DO iy=1,ny
                     ipt=(iy-1)*nx+ix
                     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
                        DO k=1,NZP1
                           VAROUT(ix,iy,k)=VEC_mean(ipt,k,i)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSE
               write(nuerr,*) 'You need to extend the outputs in'//
     &              'write_means'
            ENDIF
c            WRITE(nuout,*) 'In output_inst for varout, ivar=',ivar
c            WRITE(nuout,*) 'Attempting to write to ncid=',mean_ncid_out,
c     &           'and varid=',mean_varid(ivar)
            status=NF_PUT_VARA_REAL(
     &           mean_ncid_out,mean_varid(ivar),start,count,VAROUT)
            IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
            i=i+1
         ENDIF
      ENDDO
      
      start(3)=start(4)
      count(3)=1     
      i=1
      DO ivar=1,N_SINGOUTS
         IF (L_MEAN_SINGOUT(ivar)) THEN
            SINGOUT(:,:) = missval
            DO ix=1,nx
               DO iy=1,ny
                  ipt=(iy-1)*nx+ix
                  IF (kpp_3d_fields%L_OCEAN(ipt))
     &                 SINGOUT(ix,iy)=SCLR_mean(ipt,i)
               ENDDO
            ENDDO            
            WRITE(nuout,*) 'In output_inst for singout, ivar=',ivar
            status=NF_PUT_VARA_REAL(
     &           mean_ncid_out,mean_singid(ivar),start,count,SINGOUT)
            IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
            i=i+1
         ENDIF
      ENDDO
      start(3)=1

c     Increment counter for time dimension of NetCDF file
      nout_mean=nout_mean+1

      RETURN
      END

      SUBROUTINE mean_output(kpp_3d_fields,VEC_mean,SCLR_mean)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
! Automatically includes parameter.inc
#include <kpp_3d_type.com>
      include 'output.com'
      include 'times.com'
c      include 'ocn_paras.com'
#include <landsea.com>
#include <relax_3d.com>
      include 'couple.com'
c      include 'kprof_in.com'
      include 'ocn_advec.com'
#include <fcorr_in.com>
#include <sfcorr_in.com>

c      REAL VEC_mean(NPTS,NZP1,NVEC_MEAN),
c     &     SCLR_mean(NPTS,NSCLR_MEAN)
      REAL,intent(inout) :: VEC_mean(NPTS,NZP1,NVEC_MEAN),
     +  SCLR_mean(NPTS,NSCLR_MEAN)
         
      TYPE(kpp_3d_type) :: kpp_3d_fields
      INTEGER i,j,k,ivar,ix,iy,ipt,ipt_globe
      
      i=1
      DO ivar=1,N_VAROUTS
         IF (L_MEAN_VAROUT(ivar)) THEN
!            WRITE(6,*) 'Doing means for ivar=',ivar,', i=',i
            IF (ivar.EQ.1) THEN                                
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean)
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%U(j,k,1) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.2) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%U(j,k,2) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.3) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%X(j,k,1) / 
     +                       ndtout_mean
     &                       + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.4) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%X(j,k,2) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.5) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean)
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%buoy(j,k) /
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO               
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.6) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%wu(j,k-1,1) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO   
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.7) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%wu(j,k-1,2) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO              
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.8) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean)
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%wX(j,k-1,1) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO              
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.9) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean)
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%wX(j,k-1,2) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.10) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean)
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%wX(j,k-1,NSP1) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.11) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean)
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%wXNT(j,k-1,1) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO               
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.12) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=2,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%difm(j,k-1) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO               
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.13) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=2,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%dift(j,k-1) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO               
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.14) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=2,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%difs(j,k-1) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.15) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%rho(j,k) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.16) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%cp(j,k) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.17) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i) = kpp_3d_fields%scorr(j,k) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.18) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZ
                        VEC_mean(j,k,i) = kpp_3d_fields%Rig(j,k) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.19) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZ
                        VEC_mean(j,k,i) = kpp_3d_fields%dbloc(j,k) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.20) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZ
                        VEC_mean(j,k,i) = kpp_3d_fields%shsq(j,k) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.21) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i)=kpp_3d_fields%Tinc_fcorr(j,k) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.22) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i)=kpp_3d_fields%ocnTcorr(j,k) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ELSEIF (ivar.EQ.23) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndtout_mean) 
!$OMP& SHARED(VEC_mean,i)
!$OMP DO SCHEDULE(static)
#endif
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j)) THEN
                     DO k=1,NZP1
                        VEC_mean(j,k,i)=kpp_3d_fields%Sinc_fcorr(j,k) / 
     +                       ndtout_mean + VEC_mean(j,k,i)
                     ENDDO
                  ENDIF
               ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            ENDIF                  
            i=i+1
         ENDIF
      ENDDO
      i=1
      DO ivar=1,N_SINGOUTS         
!         WRITE(6,*) 'Means with ivar=',ivar,'i=',i
         IF (L_MEAN_SINGOUT(ivar)) THEN
!            WRITE(6,*) 'Doing means for scalar ivar=',ivar,', i=',i
            IF (ivar.EQ.1) THEN
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j))  
     &                 SCLR_mean(j,i) = kpp_3d_fields%hmix(j) / 
     +                 ndtout_mean + SCLR_mean(j,i)
               ENDDO
            ELSEIF (ivar.EQ.2) THEN
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j))
     &                 SCLR_mean(j,i) = kpp_3d_fields%fcorr(j) / 
     +                 ndtout_mean + SCLR_mean(j,i)
               ENDDO
            ELSEIF (ivar.EQ.3) THEN 
		DO j=1,NPTS
                   IF (kpp_3d_fields%L_OCEAN(j))
     &                  SCLR_mean(j,i) = kpp_3d_fields%sflux(j,1,5,0) / 
     +                  ndtout_mean + SCLR_mean(j,i)
			
		ENDDO
            ELSEIF (ivar.EQ.4) THEN
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j))
     &                 SCLR_mean(j,i) = kpp_3d_fields%sflux(j,2,5,0) / 
     +                 ndtout_mean + SCLR_mean(j,i)
               ENDDO 
	    ELSEIF (ivar.EQ.5) THEN
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j))
     &                 SCLR_mean(j,i) = kpp_3d_fields%sflux(j,3,5,0) / 
     +                 ndtout_mean + SCLR_mean(j,i)
               ENDDO
            ELSEIF (ivar.EQ.6) THEN
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j))
     &                 SCLR_mean(j,i) = kpp_3d_fields%sflux(j,4,5,0) / 
     +                 ndtout_mean + SCLR_mean(j,i)
               ENDDO
            ELSEIF (ivar .EQ. 7) THEN
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j))
     &                 SCLR_mean(j,i) = kpp_3d_fields%sflux(j,6,5,0) / 
     +                 ndtout_mean + SCLR_mean(j,i)
               ENDDO
            ELSEIF (ivar.EQ.8) THEN
               DO ix=ifirst,ilast
                  DO iy=jfirst,jlast
                     ipt_globe=(iy-1)*NX_GLOBE+ix
                     ipt=(iy-jfirst)*nx+(ix-ifirst+1)
                     SCLR_mean(ipt,i)=SCLR_mean(ipt,i) + 
     &                    kpp_3d_fields%cplwght(ipt_globe) / ndtout_mean
                  ENDDO 
               ENDDO
            ELSEIF (ivar.EQ.9) THEN
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j))
     &                 SCLR_mean(j,i)=kpp_3d_fields%freeze_flag(j) /
     +                 ndtout_mean + SCLR_mean(j,i)
               ENDDO
            ELSEIF (ivar.EQ.10) THEN
               DO j=1,NPTS
                  IF (kpp_3d_fields%L_OCEAN(j))
     &                 SCLR_mean(j,i)=kpp_3d_fields%reset_flag(j) /
     +                 ndtout_mean + SCLR_mean(j,i)
               ENDDO
            ENDIF
            i=i+1
         ENDIF
      ENDDO
            
      RETURN
      END

c$$$      SUBROUTINE mean_vec(kpp_3d_fields,VEC,VEC_mean)
c$$$      
c$$$      IMPLICIT NONE
c$$$
c$$$! Automatically includes parameter.inc!      
c$$$#include <kpp_3d_type.com>
c$$$      include 'output.com'
c$$$#include <landsea.com>
c$$$
c$$$      TYPE(kpp_3d_type) :: kpp_3d_fields
c$$$      REAL,intent(in) :: VEC(:,:)
c$$$      REAL,intent(inout) :: VEC_mean(:,:)
c$$$      INTEGER i,len_vec,len_vec_mean
c$$$      
c$$$      WRITE(6,*) 'Mean vec is computing sizes'
c$$$      len_vec=SIZE(VEC(1,:))
c$$$      len_vec_mean=SIZE(VEC_mean(1,:))      
c$$$
c$$$      IF (len_vec.ne.len_vec_mean) WRITE(6,*) 'Error: VEC and VEC_mean',
c$$$     &     'do not have the same number of vertical levels.'
c$$$      WRITE(6,*) 'len_vec=',len_vec,', len_vec_mean=',len_vec_mean
c$$$      
c$$$      DO i=1,NPTS
c$$$         IF (kpp_3d_fields%L_OCEAN(i))
c$$$     &        VEC_mean(i,:) = VEC_mean(i,:) + VEC(i,:) / ndtout_mean
c$$$      ENDDO
c$$$
c$$$      RETURN
c$$$      END

c$$$      SUBROUTINE mean_sclr(kpp_3d_fields,SCLR,SCLR_mean)
c$$$      
c$$$      IMPLICIT NONE
c$$$
c$$$! Automatically includes parameter.inc
c$$$#include <kpp_3d_type.com>
c$$$      include 'output.com'
c$$$#include <landsea.com>
c$$$
c$$$      TYPE(kpp_3d_type) :: kpp_3d_fields
c$$$      REAL SCLR(NPTS),SCLR_mean(NPTS)
c$$$      INTEGER i
c$$$      
c$$$      DO i=1,NPTS
c$$$         IF (kpp_3d_fields%L_OCEAN(i))
c$$$     &        SCLR_mean(i) = SCLR_mean(i) + SCLR(i) / ndtout_mean
c$$$
c$$$      ENDDO
c$$$         
c$$$      RETURN
c$$$      END

      SUBROUTINE output_close
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
      include 'output.com'

      INTEGER status

      status=NF_CLOSE(ncid_out)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      RETURN
      END

      SUBROUTINE output_open
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
      include 'output.com'

      INTEGER status

      status=NF_OPEN(output_file,NF_WRITE,ncid_out)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      RETURN
      END

      SUBROUTINE mean_output_open

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
      include 'netcdf.inc'
      include 'output.com'
      INTEGER status

      status=NF_OPEN(mean_output_file,NF_WRITE,mean_ncid_out)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(6,*) 'Mean output file opened with ncid=',mean_ncid_out
      
      RETURN
      END

      SUBROUTINE mean_output_close

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
      include 'netcdf.inc'
      include 'output.com'
      INTEGER status

      status=NF_CLOSE(mean_ncid_out)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      RETURN
      END

      SUBROUTINE mean_init_output(kpp_3d_fields,kpp_const_fields)

      IMPLICIT NONE
      
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
#include <kpp_3d_type.com>
      include 'output.com'
      INTEGER status
      
      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

      INTEGER temp_varid(N_VAROUTS),temp_singid(N_SINGOUTS)
      LOGICAL temp_varout(N_VAROUTS),temp_singout(N_SINGOUTS)

      nout_mean=1

      temp_varout=L_VAROUT
      L_VAROUT=L_MEAN_VAROUT
      temp_singout=L_SINGOUT
      L_SINGOUT=L_MEAN_SINGOUT
      temp_varid=varid
      temp_singid=singid

      CALL init_output(mean_output_file,mean_ncid_out,kpp_3d_fields,
     +     kpp_const_fields)

      mean_varid = varid
      mean_singid = singid
      varid = temp_varid
      singid = temp_singid
      L_VAROUT = temp_varout
      L_SINGOUT = temp_singout
            

      RETURN
      END

      
