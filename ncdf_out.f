      SUBROUTINE init_output(filename,ncid,kpp_3d_fields,
     +     kpp_const_fields,dt_vec,dt_sing,vec_varids,sing_varids,
     +     L_MEAN,L_INST)
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
c#include <landsea.com>
      include 'vert_pgrid.com'
      include 'output.com'
      include 'times.com'
      include 'initialcon.com'
c      include 'location.com'

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      INTEGER i,j,k,m,ix,iy,ncid,status,extra_time
      CHARACTER*11 varname(N_VAROUTS),singname(N_SINGOUTS)
      CHARACTER*50 longname(N_VAROUTS),singlong(N_SINGOUTS),
     &     filename
      CHARACTER*15 units(N_VAROUTS),singunits(N_SINGOUTS)
      CHARACTER*6 type
      CHARACTER*7 time_name
      REAL*4, allocatable :: time_out(:)
      INTEGER :: vec_varids(N_VAROUTS),sing_varids(N_SINGOUTS),
     &     dt_vec(N_VAROUTS),dt_sing(N_VAROUTS),n_tdims,
     &     vert_dimid(N_VAROUTS),dims(4),zflag(N_VAROUTS), ! 1 for z-levels, 0 for d-levels      
     &     time_dimids(N_VAROUTS+N_SINGOUTS), ! Contains a value for each diagnostic
     &     time_varids(N_VAROUTS+N_SINGOUTS), ! Contains only unique time variables
     &     dt_timeids(N_VAROUTS+N_SINGOUTS)
      REAL*4 ZOUT(NZP1),ALON(NX),ALAT(NY)
      REAL*4 delta
      LOGICAL tdim_found,L_MEAN,L_INST        ! If L_MEAN, time values are output as midpoints
                                              ! If L_INST, an extra value is added for initial conditions
                                              ! if .NOT. L_RESTART.                                              
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
     &     'comp_flag',
     &     'dampu_flag',
     &     'dampv_flag'/
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
     &     'Number of integrations (<0 = isothermal reset)',
     &     'Fraction of levels with ui~u**2',
     &     'Fraction of levels with vi~v**2'/

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
     &     'unitless',
     &     'fraction',
     &     'fraction'/
      
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
      delta=0.0
      IF (NY .GT. 1) delta=alat(2)-alat(1)
      CALL MY_NCDF_DEF_DIM (
     &     ncid,latdim,ny,lat_id,'latitude','deg',delta,' ')
      CALL MY_NCDF_DEF_DIM (
     &     ncid,zdim,NZP1,z_id,'z','m',0.0,' ')
      CALL MY_NCDF_DEF_DIM (
     &     ncid,hdim,NZP1,h_id,'h','m',0.0,'Layer Thickness')
      CALL MY_NCDF_DEF_DIM (
     &     ncid,ddim,NZP1,d_id,'d','m',0.0,'Depth of Interfaces')
      delta=dtout
      
      IF (dt_vec(1).gt.0) THEN
         n_tdims=1
         CALL MY_NCDF_DEF_DIM(ncid,time_dimids(1),NF_UNLIMITED,
     +        time_varids(1),'time_1','days',delta,' ')
      ELSE
         n_tdims=0
      ENDIF         
      IF (L_INST .AND. .NOT. L_RESTART) THEN
         ! Add an extra point to time dimensions for initial condition
         extra_time=1
      ELSE
         extra_time=0
      ENDIF

      DO j=2,N_VAROUTS
         tdim_found=.FALSE.
         DO k=1,j-1 
            IF (dt_vec(j) .eq. dt_vec(k)) THEN
               tdim_found=.TRUE.
               time_dimids(j)=time_dimids(k)
            ENDIF
         ENDDO
         IF (.NOT.tdim_found) THEN 
            n_tdims=n_tdims+1
            dt_timeids(n_tdims)=dt_vec(j)
            WRITE(time_name,'(A5,I0)') 'time_',n_tdims
            IF (n_tdims .eq. 1) THEN 
               CALL MY_NCDF_DEF_DIM(ncid,time_dimids(j),NF_UNLIMITED,
     +              time_varids(1),'time_1','days',delta,' ')
            ELSE
               CALL MY_NCDF_DEF_DIM (ncid,time_dimids(j),
     +              ndt_per_file/dt_vec(j)+extra_time,
     &              time_varids(n_tdims),time_name,'days',delta,' ')
            ENDIF
         ENDIF
      ENDDO

      DO j=1,N_SINGOUTS
         tdim_found=.FALSE.
         DO k=1,N_VAROUTS
            IF (dt_sing(j) .eq. dt_vec(k)) THEN
               tdim_found=.TRUE.
               time_dimids(j+N_VAROUTS)=time_dimids(k)
               WRITE(6,*) j,k,time_dimids(j+N_VAROUTS)
            ENDIF
         ENDDO
         IF (j .gt. 1 .and. .NOT. tdim_found) THEN
            DO k=1,j-1
               IF (dt_sing(j) .eq. dt_sing(k)) THEN
                  tdim_found=.TRUE.
                  time_dimids(j+N_VAROUTS)=time_dimids(k+N_VAROUTS)
               ENDIF
            ENDDO
         ENDIF
         IF (.NOT. tdim_found) THEN 
            n_tdims=n_tdims+1
            dt_timeids(n_tdims)=dt_sing(j)
            WRITE(time_name,'(A5,I0)') 'time_',n_tdims
            IF (n_tdims .eq. 1) THEN
               CALL MY_NCDF_DEF_DIM(ncid,time_dimids(j+N_VAROUTS),
     +              NF_UNLIMITED,time_varids(1),'time_1','days',
     +              delta,' ')
            ELSE
               CALL MY_NCDF_DEF_DIM(ncid,time_dimids(j+N_VAROUTS),
     +              ndt_per_file/dt_sing(j)+extra_time,
     +              time_varids(n_tdims),time_name,'days',delta,' ')
            ENDIF
         ENDIF
      ENDDO

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
      DO k=1,N_VAROUTS
         dims(3)=vert_dimid(k)
         dims(4)=time_dimids(k)
         IF (dt_vec(k) .gt. 0)
     &        CALL MY_NCDF_DEF_VAR (
     &        ncid,vec_varids(k),4,dims,varname(k),units(k),
     &        longname(k))
      ENDDO
      
      WRITE(6,*) time_dimids
      DO k=1,N_SINGOUTS
         dims(3)=time_dimids(N_VAROUTS+k)
         IF (dt_sing(k) .gt. 0)
     &        CALL MY_NCDF_DEF_VAR (
     &        ncid,sing_varids(k),3,dims,singname(k),singunits(k),
     &        singlong(k))
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

      WRITE(6,*) 'Number of time dimensions = ',n_tdims      
      DO k=1,n_tdims         
         allocate(time_out(ndt_per_file/dt_timeids(k)+extra_time))
         DO i=1-extra_time,ndt_per_file/dt_timeids(k)
c Need to write midpoints for means but output times for instaneous
            time_out(i+extra_time)=i*kpp_const_fields%dto/
     +           kpp_const_fields%spd*dt_timeids(k)
            IF (L_MEAN) 
     +           time_out(i+extra_time)=time_out(i+extra_time)-
     +           0.5*kpp_const_fields%dto/kpp_const_fields%spd*
     +           dt_timeids(k)
         ENDDO
         WRITE(6,*) 'Writing time_out=',time_out(:)
         WRITE(6,*) 'ndt_per_file=',ndt_per_file,'time_varids(k)=',
     +        time_varids(k)
         status=NF_PUT_VARA_REAL(ncid,time_varids(k),
     +        (/1/),(/ndt_per_file/dt_timeids(k)+extra_time/),
     +        time_out)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         deallocate(time_out)
      ENDDO

      write(nuout,*) 'NCDF output file initialized with name ',
     +      filename
      
      status=NF_CLOSE(ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      RETURN
      END

      SUBROUTINE output_inst(kpp_3d_fields,kpp_const_fields,diag_num,
     +     varid,ntout)
      
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

      REAL*4, allocatable :: varout(:,:,:), singout(:,:), 
     +     temp_2d(:,:),temp_1d(:)
      REAL*4 TOUT
      
      INTEGER start(4),count(4),status
      INTEGER k,ivar,diag_num,varid,ntout
      INTEGER ix,iy,ipt
 
      allocate(varout(NX,NY,NZP1))
      allocate(singout(NX,NY))
      allocate(temp_2d(NPTS,NZP1))
      allocate(temp_1d(NPTS))

      count(1)=NX
      count(2)=NY
      count(3)=NZP1
      count(4)=1
      
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=ntout

c     NPK 25/2/08.
c     Stop opening and closing output files with each flush of output.
c     Just do it once at the beginning/end of the simulation.
c      call output_open
c      TOUT=kpp_const_fields%time
c      status=NF_PUT_VAR1_REAL(ncid_out,time_id,nout,TOUT)

c      write(nuout,*) 'Writing output at timestep ',ntime+nstart,
c     +        ' Time=',TOUT

      IF (diag_num .le. N_VAROUTS) THEN 
         SELECT CASE (diag_num)
         CASE (1)
            temp_2d(:,:)=kpp_3d_fields%U(:,:,1)
         CASE (2)
            temp_2d(:,:)=kpp_3d_fields%U(:,:,2)
         CASE (3)
            temp_2d(:,:)=kpp_3d_fields%X(:,:,1)
         CASE (4)
            DO k=1,NZP1
               temp_2d(:,k)=kpp_3d_fields%X(:,k,2)+
     +              kpp_3d_fields%Sref(:)
            ENDDO
         CASE(5)
            temp_2d(:,:)=kpp_3d_fields%buoy(:,1:NZP1)
         CASE(6)
            temp_2d(:,:)=kpp_3d_fields%wU(:,0:NZ,1)
         CASE(7)
            temp_2d(:,:)=kpp_3d_fields%wU(:,0:NZ,2)
         CASE(8)
            temp_2d(:,:)=kpp_3d_fields%wX(:,0:NZ,1)
         CASE(9)
            temp_2d(:,:)=kpp_3d_fields%wX(:,0:NZ,2)
         CASE(10)
            temp_2d(:,:)=kpp_3d_fields%wX(:,0:NZ,NSP1)
         CASE(11)
            temp_2d(:,:)=kpp_3d_fields%wXNT(:,0:NZ,1)
         CASE(12)
            temp_2d(:,1)=0.0
            temp_2d(:,2:NZP1)=kpp_3d_fields%difm(:,1:NZ)
         CASE(13)
            temp_2d(:,1)=0.0
            temp_2d(:,2:NZP1)=kpp_3d_fields%dift(:,1:NZ)
         CASE(14)
            temp_2d(:,1)=0.0
            temp_2d(:,2:NZP1)=kpp_3d_fields%difs(:,1:NZ)
         CASE(15)
            temp_2d(:,:)=kpp_3d_fields%rho(:,1:NZP1)
         CASE(16)
            temp_2d(:,:)=kpp_3d_fields%cp(:,1:NZP1)
         CASE(17)
            temp_2d(:,:)=kpp_3d_fields%scorr(:,:)
         CASE(18)
            temp_2d(:,:)=kpp_3d_fields%Rig(:,:)
         CASE(19)
            temp_2d(:,1:NZ)=kpp_3d_fields%dbloc(:,1:NZ)
            temp_2d(:,NZP1)=0.0
         CASE(20)
            temp_2d(:,:)=kpp_3d_fields%Shsq(:,:)
         CASE(21)
            temp_2d(:,:)=kpp_3d_fields%Tinc_fcorr(:,:)
         CASE(22)
            temp_2d(:,:)=kpp_3d_fields%ocnTcorr(:,:)
         CASE(23)
            temp_2d(:,:)=kpp_3d_fields%Sinc_fcorr(:,:)
c         CASE DEFAULT
c            WRITE(6,*) 'You need to add more outputs in OUTPUT_INST'            
         END SELECT
      
         CALL REFORMAT_MASK_OUTPUT_2D(temp_2d,kpp_3d_fields%L_OCEAN,
     +        missval,varout)
      
         status=NF_PUT_VARA_REAL(
     &        ncid_out,varid,start,count,VAROUT)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      ELSE
         start(3)=start(4)
         count(3)=1
         SELECT CASE (diag_num-N_VAROUTS)
         CASE (1)
            temp_1d(:)=kpp_3d_fields%hmix(:)
         CASE (2)
            temp_1d(:)=kpp_3d_fields%fcorr(:)
         CASE (3)
            temp_1d(:)=kpp_3d_fields%sflux(:,1,5,0)
         CASE (4)
            temp_1d(:)=kpp_3d_fields%sflux(:,2,5,0)
         CASE (5)
            temp_1d(:)=kpp_3d_fields%sflux(:,3,5,0)
         CASE (6)
            temp_1d(:)=kpp_3d_fields%sflux(:,4,5,0)
         CASE (7)
            temp_1d(:)=kpp_3d_fields%sflux(:,6,5,0)
         CASE (8)
            DO ix=ifirst,ilast
               DO iy=jfirst,jlast
                  ipt=(iy-1)*NX_GLOBE+ix
                  temp_1d((ix-ifirst)*NY+iy-jfirst+1)=
     +                 kpp_3d_fields%cplwght(ipt)                     
               ENDDO
            ENDDO
         CASE (9)
            temp_1d(:)=kpp_3d_fields%freeze_flag(:)
         CASE (10)
            temp_1d(:)=kpp_3d_fields%reset_flag(:)
         CASE (11)
            temp_1d(:)=kpp_3d_fields%dampu_flag(:)
         CASE (12)
            temp_1d(:)=kpp_3d_fields%dampv_flag(:)
         CASE DEFAULT
            WRITE(6,*) 'You need to add more outputs in '//
     +           'OUTPUT_INST'
         END SELECT
            
         CALL REFORMAT_MASK_OUTPUT_1D(temp_1d,kpp_3d_fields%L_OCEAN,
     +        missval,singout)
         
         status=NF_PUT_VARA_REAL(
     &        ncid_out,varid,start,count,SINGOUT)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      ENDIF
      start(3)=1
      
      ntout=ntout+1
c     NPK 25/2/08
c     Stop opening and closing output files with each flush of output.
c     Just do it once at the beginning/end of the simulation.
c      call output_close
      WRITE(nuout,*) ' Output successfully written'
      RETURN
      END

      SUBROUTINE reformat_mask_output_1d(oned_in,mask,missval,
     +     twod_out)
      IMPLICIT NONE
#include <parameter.inc>

      REAL*4,intent(in) :: oned_in(NPTS),missval
      REAL*4,intent(out) :: twod_out(NX,NY)
      LOGICAL,intent(in) :: mask(NPTS)
      INTEGER :: i,j,ipt

      DO i=1,NX
         DO j=1,NY
            ipt=(j-1)*NX+i
            IF (mask(ipt)) THEN
               twod_out(i,j)=oned_in(ipt)
            ELSE
               twod_out(i,j)=missval
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE reformat_mask_output_1d

      SUBROUTINE reformat_mask_output_2d(twod_in,mask,missval,
     +     threed_out)

      IMPLICIT NONE
#include <parameter.inc>

      REAL*4,intent(in) :: twod_in(NPTS,NZP1),missval
      REAL*4,intent(out) :: threed_out(NX,NY,NZP1)
      LOGICAL,intent(in) :: mask(NPTS)
      INTEGER :: i,j,ipt

      DO i=1,NX
         DO j=1,NY
            ipt=(j-1)*NX+i
            IF (mask(ipt)) THEN
               threed_out(i,j,:)=twod_in(ipt,:)
            ELSE
               threed_out(i,j,:)=missval
            ENDIF
         ENDDO
      ENDDO
                  
      RETURN
      END SUBROUTINE reformat_mask_output_2d

      SUBROUTINE write_means(kpp_3d_fields,kpp_const_fields,
     +     VEC_mean,SCLR_mean,diag_num,varid,ntout)
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
      REAL*4,allocatable :: VAROUT(:,:,:), SINGOUT(:,:),
     +     temp_2d(:,:)
      REAL*4 TOUT

      INTEGER i,ivar,ipt,ix,iy,start(4),count(4),k,status,
     +     diag_num,varid,ntout,mean_num

      allocate(VAROUT(NX,NY,NZP1))
      allocate(SINGOUT(NX,NY))
      allocate(temp_2d(NPTS,NZP1))
c      TOUT=kpp_const_fields%time-(ndtout_mean*kpp_const_fields%dtsec/
c     +     (ndtocn*kpp_const_fields%spd))*0.5
c      status=NF_PUT_VAR1_REAL(mean_ncid_out,time_id,nout_mean,TOUT)
c      IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
c      write(nuout,*) 'Writing mean output at timestep ',ntime+nstart,
c     +     ' Time=',TOUT,' time_id = ',time_id
      
      count(1)=NX
      count(2)=NY
      count(3)=NZP1
      count(4)=1
      
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=ntout
      
      IF (diag_num .le. N_VAROUTS) THEN 
         mean_num=1
         DO i=1,diag_num-1
            IF (ndt_varout_mean(i) .gt. 0) 
     +           mean_num=mean_num+1
         ENDDO
         SELECT CASE (diag_num)
         CASE (4) 
            DO k=1,NZP1
               temp_2d(:,k)=VEC_mean(:,k,mean_num)+
     +              kpp_3d_fields%Sref(:)
            ENDDO
         CASE (12,13,14)
            temp_2d(:,1)=0.
            temp_2d(:,2:NZP1)=VEC_mean(:,2:NZP1,mean_num)
         CASE (18,19)
            temp_2d(:,1:NZ)=VEC_mean(:,1:NZ,mean_num)
            temp_2d(:,NZP1)=0.               
         CASE DEFAULT
            temp_2d(:,:)=VEC_mean(:,:,mean_num)
         END SELECT
         WRITE(6,*) 'Calling reformat_mask_output_2d for ivar=',ivar
         CALL reformat_mask_output_2d(temp_2d,kpp_3d_fields%L_OCEAN,
     +        missval,varout)
         
         WRITE(6,*) 'Writing to ncid=',mean_ncid_out,' varid=',
     +        varid,' with start =',start,' and count =',
     +        count
         status=NF_PUT_VARA_REAL(
     &        mean_ncid_out,varid,start,count,VAROUT)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         VEC_mean(:,:,mean_num)=0.
         i=i+1         
      ELSE
         start(3)=start(4)
         count(3)=1     
         mean_num=1
         DO i=1,(diag_num-N_VAROUTS-1)
            IF (ndt_singout_mean(i) .gt. 0)
     +           mean_num=mean_num+1
         ENDDO
         SINGOUT(:,:) = missval
         DO ix=1,nx
            DO iy=1,ny
               ipt=(iy-1)*nx+ix
               IF (kpp_3d_fields%L_OCEAN(ipt))
     &              SINGOUT(ix,iy)=SCLR_mean(ipt,mean_num)
            ENDDO
         ENDDO            
         WRITE(nuout,*) 'In write_means for singout, ivar=',ivar
         status=NF_PUT_VARA_REAL(
     &        mean_ncid_out,varid,start,count,SINGOUT)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         SCLR_mean(:,mean_num)=0.
         i=i+1
      ENDIF     
      start(3)=1

c     Increment counter for time dimension of NetCDF file
      ntout=ntout+1

      RETURN
      END

      SUBROUTINE output_range(kpp_3d_fields,kpp_const_fields,
     +     VEC_range,SCLR_range,diag_num,varid,ntout)
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

      REAL :: VEC_range(NPTS,NZP1,NVEC_RANGE,2),
     +     SCLR_range(NPTS,NSCLR_RANGE,2)
      REAL*4,allocatable :: VAROUT(:,:,:), SINGOUT(:,:),
     +     temp_2d(:,:)
      REAL*4 TOUT

      INTEGER i,j,ivar,ipt,ix,iy,start(4),count(4),k,status,
     +     diag_num,varid,ntout,range_num

      allocate(VAROUT(NX,NY,NZP1))
      allocate(SINGOUT(NX,NY))
      allocate(temp_2d(NPTS,NZP1))
      
      count(1)=NX
      count(2)=NY
      count(3)=NZP1
      count(4)=1
      
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=ntout
      
      DO j=1,2
         WRITE(6,*) 'j=',j,'diag_num=',diag_num
         IF (diag_num .le. N_VAROUTS) THEN 
            range_num=1
            DO i=1,diag_num-1
               IF (ndt_varout_range(i) .gt. 0) 
     +              range_num=range_num+1
            ENDDO
            WRITE(6,*) 'range_num=',range_num
            SELECT CASE (diag_num)
            CASE (4) 
               DO k=1,NZP1
                  temp_2d(:,k)=VEC_range(:,k,range_num,j)+
     +                 kpp_3d_fields%Sref(:)
               ENDDO
            CASE (12,13,14)
               temp_2d(:,1)=0.
               temp_2d(:,2:NZP1)=VEC_range(:,2:NZP1,range_num,j)
            CASE (18,19)
               temp_2d(:,1:NZ)=VEC_range(:,1:NZ,range_num,j)
               temp_2d(:,NZP1)=0.               
            CASE DEFAULT
               temp_2d(:,:)=VEC_range(:,:,range_num,j)
            END SELECT
            WRITE(6,*) 'Calling reformat_mask_output_2d diag=',diag_num
            CALL reformat_mask_output_2d(temp_2d,kpp_3d_fields%L_OCEAN,
     +           missval,varout)           
            SELECT CASE (j)
            CASE (1)                              
               WRITE(6,*) 'Writing to ncid=',min_ncid_out,' varid=',
     +              varid,' with start =',start,' and count =',
     +              count
               status=NF_PUT_VARA_REAL(
     &              min_ncid_out,varid,start,count,VAROUT)
               IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
               VEC_range(:,:,range_num,j)=2e20
            CASE (2)
               WRITE(6,*) 'Writing to ncid=',max_ncid_out,' varid=',
     +              varid,' with start =',start,' and count =',
     +              count
               status=NF_PUT_VARA_REAL(
     &              max_ncid_out,varid,start,count,VAROUT)
               IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
               VEC_range(:,:,range_num,j)=-2e20
            END SELECT                     
c     i=i+1         
         ELSE
            start(3)=start(4)
            count(3)=1     
            range_num=1
            DO i=1,(diag_num-N_VAROUTS-1)
               IF (ndt_singout_range(i) .gt. 0)
     +              range_num=range_num+1
            ENDDO
            SINGOUT(:,:) = missval
            DO ix=1,nx
               DO iy=1,ny
                  ipt=(iy-1)*nx+ix
                  IF (kpp_3d_fields%L_OCEAN(ipt))
     &                 SINGOUT(ix,iy)=SCLR_range(ipt,range_num,j)
               ENDDO
            ENDDO    
            SELECT CASE (j)
            CASE (1)               
               status=NF_PUT_VARA_REAL(
     &              min_ncid_out,varid,start,count,SINGOUT)
               IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
               SCLR_range(:,range_num,j)=2e20
            CASE (2)
               status=NF_PUT_VARA_REAL(
     &              max_ncid_out,varid,start,count,SINGOUT)
               IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
               SCLR_range(:,range_num,j)=-2e20
            END SELECT
c     i=i+1
         ENDIF     
      ENDDO
      start(3)=1

c     Increment counter for time dimension of NetCDF file
      ntout=ntout+1

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
#include <landsea.com>
#include <relax_3d.com>
      include 'couple.com'
      include 'ocn_advec.com'
#include <fcorr_in.com>
#include <sfcorr_in.com>

      REAL,intent(inout) :: VEC_mean(NPTS,NZP1,NVEC_MEAN),
     +  SCLR_mean(NPTS,NSCLR_MEAN)
      REAL, allocatable :: field(:,:),vec(:)

      TYPE(kpp_3d_type) :: kpp_3d_fields
      INTEGER i,j,k,ivar,ix,iy,ipt,ipt_globe,upper_limit,lower_limit
      
      allocate(field(NPTS,NZP1))
      allocate(vec(NPTS))

      i=1
      DO ivar=1,N_VAROUTS
         IF (ndt_varout_mean(ivar) .gt. 0) THEN
            SELECT CASE (ivar)
            CASE(1)
               field(:,:)=kpp_3d_fields%U(:,:,1)
            CASE(2)
               field(:,:)=kpp_3d_fields%U(:,:,2)
            CASE(3)
               field(:,:)=kpp_3d_fields%X(:,:,1)
            CASE(4)
               field(:,:)=kpp_3d_fields%X(:,:,2)
            CASE(5)
               field(:,:)=kpp_3d_fields%buoy(:,1:NZP1)
            CASE(6)
               field(:,:)=kpp_3d_fields%wu(:,0:NZ,1)
            CASE(7)
               field(:,:)=kpp_3d_fields%wu(:,0:NZ,2)
            CASE(8)
               field(:,:)=kpp_3d_fields%wx(:,0:NZ,1)
            CASE(9)
               field(:,:)=kpp_3d_fields%wx(:,0:NZ,2)
            CASE(10)
               field(:,:)=kpp_3d_fields%wx(:,0:NZ,NSP1)
            CASE(11)
               field(:,:)=kpp_3d_fields%wXNT(:,0:NZ,1)
            CASE(12)
               field(:,:)=kpp_3d_fields%difm(:,1:NZP1)
            CASE(13)
               field(:,:)=kpp_3d_fields%dift(:,1:NZP1)
            CASE(14)
               field(:,:)=kpp_3d_fields%difs(:,1:NZP1)
            CASE(15)
               field(:,:)=kpp_3d_fields%rho(:,1:NZP1)
            CASE(16)
               field(:,:)=kpp_3d_fields%cp(:,1:NZP1)
            CASE(17)
               field(:,:)=kpp_3d_fields%scorr(:,:)
            CASE(18)
               field(:,:)=kpp_3d_fields%Rig(:,:)
            CASE(19)
               field(:,1:NZ)=kpp_3d_fields%dbloc(:,1:NZ)
               field(:,NZP1)=0.
            CASE(20)
               field(:,:)=kpp_3d_fields%Shsq(:,:)
            CASE(21)
               field(:,:)=kpp_3d_fields%Tinc_fcorr(:,:)
            CASE(22)
               field(:,:)=kpp_3d_fields%ocnTcorr(:,:)
            CASE(23)
               field(:,:)=kpp_3d_fields%sinc_fcorr(:,:)
            END SELECT
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndt_varout_mean)
!$OMP& SHARED(VEC_mean,i,field,upper_limit,lower_limit)
!$OMP DO SCHEDULE(static)
#endif
            DO j=1,NPTS
               IF (kpp_3d_fields%L_OCEAN(j)) THEN
                  DO k=1,NZP1
                     VEC_mean(j,k,i)=field(j,k) / ndt_varout_mean(ivar) 
     +                    + VEC_mean(j,k,i)
                  ENDDO
               ENDIF
            ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif         
            i=i+1
         ENDIF
      ENDDO
      i=1
      DO ivar=1,N_SINGOUTS         
!         WRITE(6,*) 'Means with ivar=',ivar,'i=',i
         IF (ndt_singout_mean(ivar) .gt. 0) THEN
            SELECT CASE (ivar)
            CASE(1)
               vec(:)=kpp_3d_fields%hmix(:)
            CASE(2)
               vec(:)=kpp_3d_fields%fcorr(:)
            CASE(3)
               vec(:)=kpp_3d_fields%sflux(:,1,5,0)
            CASE(4)
               vec(:)=kpp_3d_fields%sflux(:,2,5,0)
            CASE(5)
               vec(:)=kpp_3d_fields%sflux(:,3,5,0)
            CASE(6)
               vec(:)=kpp_3d_fields%sflux(:,4,5,0)
            CASE(7)
               vec(:)=kpp_3d_fields%sflux(:,6,5,0)
            CASE(8)
               DO ix=ifirst,ilast
                  DO iy=jfirst,jlast
                     ipt_globe=(iy-1)*NX_GLOBE+ix
                     ipt=(iy-jfirst)*nx+(ix-ifirst+1)
                     vec(ipt)=kpp_3d_fields%cplwght(ipt_globe)
                  ENDDO
               ENDDO
            CASE(9)
               vec(:)=kpp_3d_fields%freeze_flag(:)
            CASE(10)
               vec(:)=kpp_3d_fields%reset_flag(:)
            CASE(11)
               vec(:)=kpp_3d_fields%dampu_flag(:)
            CASE(12)
               vec(:)=kpp_3d_fields%dampv_flag(:)
            END SELECT
            DO j=1,NPTS
               IF (kpp_3d_fields%L_OCEAN(j)) 
     +              SCLR_mean(j,i)=vec(j)/ndt_singout_mean(ivar)+
     +              SCLR_mean(j,i)
            ENDDO
         ENDIF
         i=i+1      
      ENDDO
      deallocate(field)
      deallocate(vec)
            
      RETURN
      END

      SUBROUTINE range_output(kpp_3d_fields,VEC_range,SCLR_range)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
! Automatically includes parameter.inc
#include <kpp_3d_type.com>
      include 'output.com'
      include 'times.com'
#include <landsea.com>
#include <relax_3d.com>
      include 'couple.com'
      include 'ocn_advec.com'
#include <fcorr_in.com>
#include <sfcorr_in.com>

      REAL,intent(inout) :: VEC_range(NPTS,NZP1,NVEC_RANGE,2),
     +  SCLR_range(NPTS,NSCLR_RANGE,2)
      REAL, allocatable :: field(:,:),vec(:)

      TYPE(kpp_3d_type) :: kpp_3d_fields
      INTEGER i,j,k,ivar,ix,iy,ipt,ipt_globe
      
      allocate(field(NPTS,NZP1))
      allocate(vec(NPTS))

      i=1
      DO ivar=1,N_VAROUTS
         IF (ndt_varout_range(ivar) .gt. 0) THEN
            SELECT CASE (ivar)
            CASE(1)
               field(:,:)=kpp_3d_fields%U(:,:,1)
            CASE(2)
               field(:,:)=kpp_3d_fields%U(:,:,2)
            CASE(3)
               field(:,:)=kpp_3d_fields%X(:,:,1)
            CASE(4)
               field(:,:)=kpp_3d_fields%X(:,:,2)
            CASE(5)
               field(:,:)=kpp_3d_fields%buoy(:,1:NZP1)
            CASE(6)
               field(:,:)=kpp_3d_fields%wu(:,0:NZ,1)
            CASE(7)
               field(:,:)=kpp_3d_fields%wu(:,0:NZ,2)
            CASE(8)
               field(:,:)=kpp_3d_fields%wx(:,0:NZ,1)
            CASE(9)
               field(:,:)=kpp_3d_fields%wx(:,0:NZ,2)
            CASE(10)
               field(:,:)=kpp_3d_fields%wx(:,0:NZ,NSP1)
            CASE(11)
               field(:,:)=kpp_3d_fields%wXNT(:,0:NZ,1)
            CASE(12)
               field(:,:)=kpp_3d_fields%difm(:,1:NZP1)
            CASE(13)
               field(:,:)=kpp_3d_fields%dift(:,1:NZP1)
            CASE(14)
               field(:,:)=kpp_3d_fields%difs(:,1:NZP1)
            CASE(15)
               field(:,:)=kpp_3d_fields%rho(:,1:NZP1)
            CASE(16)
               field(:,:)=kpp_3d_fields%cp(:,1:NZP1)
            CASE(17)
               field(:,:)=kpp_3d_fields%scorr(:,:)
            CASE(18)
               field(:,:)=kpp_3d_fields%Rig(:,:)
            CASE(19)
               field(:,1:NZ)=kpp_3d_fields%dbloc(:,1:NZ)
               field(:,NZP1)=0.
            CASE(20)
               field(:,:)=kpp_3d_fields%Shsq(:,:)
            CASE(21)
               field(:,:)=kpp_3d_fields%Tinc_fcorr(:,:)
            CASE(22)
               field(:,:)=kpp_3d_fields%ocnTcorr(:,:)
            CASE(23)
               field(:,:)=kpp_3d_fields%sinc_fcorr(:,:)
            END SELECT
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields)
!$OMP& SHARED(VEC_range,i,field)
!$OMP DO SCHEDULE(static)
#endif
            DO j=1,NPTS
               IF (kpp_3d_fields%L_OCEAN(j)) THEN
                  DO k=1,NZP1
                     IF (field(j,k) .lt. VEC_range(j,k,i,1))
     +                    VEC_range(j,k,i,1)=field(j,k)
                     IF (field(j,k) .gt. VEC_range(j,k,i,2))
     +                    VEC_range(j,k,i,2)=field(j,k)
                  ENDDO
               ENDIF
            ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif         
            i=i+1
         ENDIF
      ENDDO
      i=1
      DO ivar=1,N_SINGOUTS         
!         WRITE(6,*) 'Means with ivar=',ivar,'i=',i
         IF (ndt_singout_mean(ivar) .gt. 0) THEN
            SELECT CASE (ivar)
            CASE(1)
               vec(:)=kpp_3d_fields%hmix(:)
            CASE(2)
               vec(:)=kpp_3d_fields%fcorr(:)
            CASE(3)
               vec(:)=kpp_3d_fields%sflux(:,1,5,0)
            CASE(4)
               vec(:)=kpp_3d_fields%sflux(:,2,5,0)
            CASE(5)
               vec(:)=kpp_3d_fields%sflux(:,3,5,0)
            CASE(6)
               vec(:)=kpp_3d_fields%sflux(:,4,5,0)
            CASE(7)
               vec(:)=kpp_3d_fields%sflux(:,6,5,0)
            CASE(8)
               DO ix=ifirst,ilast
                  DO iy=jfirst,jlast
                     ipt_globe=(iy-1)*NX_GLOBE+ix
                     ipt=(iy-jfirst)*nx+(ix-ifirst+1)
                     vec(ipt)=kpp_3d_fields%cplwght(ipt_globe)
                  ENDDO
               ENDDO
            CASE(9)
               vec(:)=kpp_3d_fields%freeze_flag(:)
            CASE(10)
               vec(:)=kpp_3d_fields%reset_flag(:)
            CASE(11)
               vec(:)=kpp_3d_fields%dampu_flag(:)
            CASE(12)
               vec(:)=kpp_3d_fields%dampv_flag(:)
            END SELECT
            DO j=1,NPTS
               IF (kpp_3d_fields%L_OCEAN(j)) THEN
                  IF (vec(j) .lt. SCLR_range(j,i,1))
     +                 SCLR_range(j,i,1)=vec(j)
                  IF (vec(j) .gt. SCLR_range(j,i,2))
     +                 SCLR_range(j,i,2)=vec(j)
               ENDIF                    
            ENDDO
         ENDIF
         i=i+1      
      ENDDO
      deallocate(field)
      deallocate(vec)
            
      RETURN
      END

      SUBROUTINE output_close(ncid)
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
      include 'output.com'

      INTEGER status
      INTEGER, intent(in) :: ncid

      status=NF_CLOSE(ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      RETURN
      END

      SUBROUTINE output_open(file,ncid)
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'netcdf.inc'
      include 'output.com'

      INTEGER status
      INTEGER,intent(out) :: ncid
      CHARACTER*50,intent(in) :: file      

      status=NF_OPEN(file,NF_WRITE,ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      RETURN
      END
      
