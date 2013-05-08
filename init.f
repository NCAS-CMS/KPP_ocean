      SUBROUTINE init_env (lstretchgrid,dscale,kpp_3d_fields,
     +     kpp_const_fields)
c     ===================
c     My init_env routine modified from ...
c     Modified  3 March 1991   -  WGL in KPP/input.f

c     the physical environment for the model:
c     vertical grid and geographic location
c
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc
#include <kpp_3d_type.com>
#include <netcdf.inc>
      include 'constants.com'
      include 'vert_pgrid.com'
      include 'location.com'
      include 'ocn_state.com'
c     Inputs
      LOGICAL lstretchgrid
      REAL dscale
      TYPE(kpp_3d_type)    :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
c     Local Variables
      REAL sumh,hsum,dfac,sk
      REAL*4 vgrid_in(NZ)
      INTEGER i,ipt,ncid,status,dimid,varid
c define vertical grid fields

      IF (L_VGRID_FILE) THEN 
         WRITE(6,*) 'Reading vertical grid from file ',vgrid_file
         status=NF_OPEN(vgrid_file,0,ncid)
         IF (status.ne.NF_NOERR) CALL HANDLE_ERR(status)
         status=NF_INQ_VARID(ncid,'d',varid)
         IF (status.ne.NF_NOERR) CALL HANDLE_ERR(status)
         status=NF_GET_VAR_REAL(ncid,varid,vgrid_in)
         IF (status.ne.NF_NOERR) CALL HANDLE_ERR(status)
         kpp_const_fields%dm(1:NZ)=vgrid_in
         status=NF_OPEN(vgrid_file,0,ncid)
         IF (status.ne.NF_NOERR) CALL HANDLE_ERR(status)
         status=NF_INQ_VARID(ncid,'h',varid)
         IF (status.ne.NF_NOERR) CALL HANDLE_ERR(status)
         status=NF_GET_VAR_REAL(ncid,varid,vgrid_in)
         kpp_const_fields%hm(1:NZ)=vgrid_in
         IF (status.ne.NF_NOERR) CALL HANDLE_ERR(status)
         status=NF_OPEN(vgrid_file,0,ncid)
         IF (status.ne.NF_NOERR) CALL HANDLE_ERR(status)
         status=NF_INQ_VARID(ncid,'z',varid)
         IF (status.ne.NF_NOERR) CALL HANDLE_ERR(status)
         status=NF_GET_VAR_REAL(ncid,varid,vgrid_in)
         kpp_const_fields%zm(1:NZ)=vgrid_in
         IF (status.ne.NF_NOERR) CALL HANDLE_ERR(status)
         status=NF_CLOSE(ncid)
         DMAX=-1.*(kpp_const_fields%zm(NZ)-kpp_const_fields%hm(NZ))
      ELSE      
         IF (lstretchgrid) THEN
            sumh = 0.0
            dfac = 1.0 - exp(-dscale)
            do 5 i = 1,(NZ)
               sk = - (float(i)-0.5)/float(NZ)
               kpp_const_fields%hm(i) = 
     +              DMAX*dfac/float(NZ)/dscale / ( 1.0 + sk*dfac )
               sumh = sumh + kpp_const_fields%hm(i)
 5          continue
         ENDIF
c     
c     layer thickness h, layer grids zgrid, interface depths d
c     
         hsum = 0.0
         kpp_const_fields%dm(0) = 0.0
         do 15 i=1,NZ
            if(lstretchgrid) then
               kpp_const_fields%hm(i) = kpp_const_fields%hm(i) 
     +              *DMAX / sumh 
            else   
               kpp_const_fields%hm(i) = DMAX / real(NZ) 
            endif   
            kpp_const_fields%zm(i) =  0.0 - (hsum + 0.5 * 
     +           kpp_const_fields%hm(i) )
            hsum = hsum + kpp_const_fields%hm(i)
            kpp_const_fields%dm(i) = hsum
c     write(nuout,*) 'h=',hm(i),' dm=',dm(i),' z=',zm(i),' i=',i
 15      continue
      ENDIF
      kpp_const_fields%hm(nzp1) = 1.e-10 
      kpp_const_fields%zm(nzp1) = -DMAX
c     
c     set fraction of land, and initial open water fraction
c     
         DO ipt=1,NPTS
            focn (ipt) = 1.
c     
c     compute geographic location
c     
c     rlat(ipt)=dlat(ipt)*twopi/360.
c     rlon(ipt)=dlon(ipt)*twopi/360.
c     
c     Coriolis Parameter ! check the necessity of this
c     
            if(abs(kpp_3d_fields%dlat(ipt)).lt.2.5) then
               kpp_3d_fields%f(ipt) = 2. * (twopi/86164.) * 
     +              sin(2.5*twopi/360.)
            else  
               kpp_3d_fields%f(ipt) = 2. * (twopi/86164.) * 
     +              sin(kpp_3d_fields%dlat(ipt)*twopi/360.)
            endif  
         ENDDO


      return
      end
*************************************************************
      subroutine init_flds (kpp_3d_fields,kpp_const_fields)
c     ====================
c     provides initial conditions 
c
c     Modified  3 March  1991  --  WGL
c
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
      include 'initialcon.com'
      include 'ocn_paras.com'
      include 'ocn_state.com'
c output
c      real X(NPTS,NZP1,NSCLR),U(NPTS,NZP1,NVEL)
      real var_in(NZP1)
c local
      INTEGER n,i,ipt

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields


      IF ( .NOT. L_INITDATA ) THEN
         write(nuerr,*) "No code for L_INITDATA=.FALSE."
           CALL MIXED_ABORT
c         DO ipt=1,npts
c            do  n=1,NVEL
c               call vprofile(var_in,SSU(ipt,n),hu(ipt),mu(ipt))
c               DO i=1,NZP1
c                  U(ipt,i,n)=var_in(i)
c               ENDDO
c            enddo
c            do n=1,NSCLR
c               call vprofile(var_in,SSX(ipt,n),hx(ipt),mx(ipt)) 
c               DO i=1,NZP1
c                  X(ipt,i,n)=var_in(i)
c               ENDDO
c            enddo 
c         ENDDO
c     
      ELSE
         call read_init(kpp_3d_fields,kpp_const_fields)
      ENDIF
c
c Calculate and remove reference salinity      
c
      DO ipt=1,npts
         kpp_3d_fields%Sref(ipt)=(kpp_3d_fields%X(ipt,1,2)+
     +        kpp_3d_fields%X(ipt,nzp1,2))/2.
         kpp_3d_fields%Ssref(ipt)=kpp_3d_fields%Sref(ipt)

c      write(nuout,1010) Sref(ipt)
c 1010 format(/,'salinity reference value =',f11.6,/)
         do i=1,nzp1
            kpp_3d_fields%X(ipt,i,2)=kpp_3d_fields%X(ipt,i,2)-
     +           kpp_3d_fields%Sref(ipt)
         enddo
c     Initial surface temp
         kpp_3d_fields%Tref(ipt) = kpp_3d_fields%X(ipt,1,1)
         IF (L_SSref) THEN
            kpp_3d_fields%Ssurf(ipt)=kpp_3d_fields%SSref(ipt)
         ELSE
            kpp_3d_fields%Ssurf(ipt)=kpp_3d_fields%X(ipt,1,2)+
     +           kpp_3d_fields%Sref(ipt)
         ENDIF
      ENDDO
      
      return
      end

************************************************************************
      SUBROUTINE vprofile(y,y0,hmixd,modeprof)
c     =================     
c     returns a vertical profile of property y via MODE
c
c     Written  3  March  1991  - WGL
c
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <parameter.inc>
      include 'constants.com'
      include 'vert_pgrid.com'
c
c Inputs
      REAL y0, y(nzp1),hmixd
      INTEGER modeprof
c Local variables
      INTEGER i
      REAL dep,alf,a,b

      if(hmixd.ge.DMAX) hmixd =  0. 99 * DMAX  
c
c STEP PROFILE
c
      if (modeprof.eq.1) then
         do 15 i=1,nz
            dep = -zm(i)
            if(dep.le.hmixd) then
               y(i) = y0
            else
               y(i) = y(nzp1)
            endif
 15      continue
c     
c     EXPONENTIAL DECREASE
c     
      else if (modeprof.eq.2) then
         alf = alog( y0 / y(nzp1) ) / ( DMAX -  hmixd )
         do 25 i=1,nz
            dep = -zm(i)
            if(dep.le.hmixd) then
               y(i) = y0
            else
               y(i) = y0 * exp( alf * (hmixd - dep) )
            endif
 25      continue
c     
c     LINEAR DECREASE
c     
      else if (modeprof.eq.3) then
         alf = (y0 - y(nzp1)) / ( hmixd - DMAX )
c     
         do 45 i = 1,nz
            dep = -zm(i)
            if(dep.le.hmixd) then
               y(i) = y0 
            else
               y(i) = y0 + alf * ( dep - hmixd )
            endif
 45      continue
c     
c     SINUSOIDAL DECREASE
c     
      else if (modeprof.eq.4) then 
         alf = ONEPI / ( DMAX - hmixd )
         a   = 0.5 * (y0 - y(nzp1) )
         b   = 0.5 * (y0 + y(nzp1) )
c     
         do 35 i = 1,nz
            dep = -zm(i)
            if(dep.le.hmixd) then
               y(i) = y0
            else
               y(i) = b + a * cos( alf * (dep - hmixd) )
            endif 
 35      continue
c     
      else
         write(nuerr,*) 'STOP in vprofile (input.f):'
         write(nuerr,*) '     unable to initialize vprofile, modeprof=',
     $        modeprof
         CALL MIXED_ABORT
      endif
c     
      RETURN
      END 
      
************************************************************************

      SUBROUTINE init_flx(kpp_3d_fields)

c     Set up parameters for calculating fluxes and initialize fluxes.
c     intermediate values computed every ndtld
c     common deltax(NJDT), xbar, denom

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
      include 'flx_profs.com'
      include 'flx_sfc.com'
      TYPE(kpp_3d_type) :: kpp_3d_fields
      integer i,ipt
c
c extrapolation parameters for fluxes
c
c
c     Initialize flux arrays
c
      call aset(kpp_3d_fields%wU,NZP1tmax*NVP1,0.0)
      call aset(kpp_3d_fields%wX,NZP1tmax*NSP1,0.0)
      call aset(kpp_3d_fields%wU,NZP1tmax*NSCLR,0.0)
      call aset(kpp_3d_fields%sflux,nsflxs*5*(njdt+1),0.0)

      DO ipt=1,npts
         do i=1,nsflxs
            kpp_3d_fields%sflux(ipt,i,5,0)=1e-20
         enddo
      ENDDO
      return
      end

      SUBROUTINE init_paras(kpp_3d_fields)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <kpp_3d_type.com>
      include 'proc_pars.com'
      
      TYPE(kpp_3d_type) :: kpp_3d_fields
      integer jerlov(npts)
      integer ipt

      LOGICAL L_PARAS

      L_PARAS=L_JERLOV

      IF (L_PARAS) THEN
         call init_parasfile
      ENDIF

      IF (L_JERLOV) THEN
         call read_ipar(kpp_3d_fields,
     &        ncid_paras,'jerlov',1,1,kpp_3d_fields%jerlov)
      ELSE
         DO ipt=1,npts
            kpp_3d_fields%jerlov(ipt)=3
         ENDDO
      ENDIF

      IF (L_PARAS) THEN
         CALL close_parasfile
      ENDIF

      RETURN
      END


 
      SUBROUTINE init_advect(kpp_3d_fields)
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
      include 'ocn_advec.com'

      TYPE(kpp_3d_type) :: kpp_3d_fields
      INTEGER nmode(npts),mode(npts,maxmodeadv)
      REAL adv(npts,maxmodeadv)

      INTEGER ipt,ivar,imode

      IF (L_ADVECT) THEN
         call init_advectfile

         call read_ipar(kpp_3d_fields,
     &        ncid_advec,'nmode_tadv',1,1,nmode)
         DO ipt=1,npts
            kpp_3d_fields%nmodeadv(ipt,1)=nmode(ipt)
         ENDDO
         call read_ipar(kpp_3d_fields,
     &        ncid_advec,'mode_tadv',maxmodeadv,1,mode)
         DO ipt=1,npts
            DO imode=1,maxmodeadv
               kpp_3d_fields%modeadv(ipt,imode,1)=mode(ipt,imode)
            ENDDO
         ENDDO
         call read_par(kpp_3d_fields,
     &        ncid_advec,'tadv',maxmodeadv,1,adv)
         DO ipt=1,npts
            DO imode=1,maxmodeadv
               kpp_3d_fields%advection(ipt,imode,1)=adv(ipt,imode)
            ENDDO
         ENDDO

         call read_ipar(kpp_3d_fields,
     &        ncid_advec,'nmode_sadv',1,1,nmode)
         DO ipt=1,npts
            kpp_3d_fields%nmodeadv(ipt,2)=nmode(ipt)
         ENDDO
         call read_ipar(kpp_3d_fields,
     &        ncid_advec,'mode_sadv',maxmodeadv,1,mode)
         DO ipt=1,npts
            DO imode=1,maxmodeadv
               kpp_3d_fields%modeadv(ipt,imode,2)=mode(ipt,imode)
            ENDDO
         ENDDO
         call read_par(kpp_3d_fields,
     &        ncid_advec,'sadv',maxmodeadv,1,adv)
         DO ipt=1,npts
            DO imode=1,maxmodeadv
               kpp_3d_fields%advection(ipt,imode,2)=adv(ipt,imode)
            ENDDO
         ENDDO

         CALL close_advectfile

      ELSE
         DO ipt=1,npts
            DO ivar=1,2
               kpp_3d_fields%nmodeadv(ipt,ivar)=0
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END

************************************************************
      SUBROUTINE init_landsea(kpp_3d_fields)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc
#include <kpp_3d_type.com>
#include <landsea.com>

      REAL landsea(npts)
      TYPE (kpp_3d_type) :: kpp_3d_fields
      
      INTEGER ipt

      IF (L_LANDSEA) THEN
         call init_landseafile

         WRITE(6,*) 'Initialized landsea file',ncid_landsea

         call read_par(kpp_3d_fields,
     &        ncid_landsea,'lsm',1,1,landsea)
         DO ipt=1,npts
            IF (landsea(ipt) .GE. 0.5) THEN
               kpp_3d_fields%L_OCEAN(ipt)=.FALSE.
            ELSE
               kpp_3d_fields%L_OCEAN(ipt)=.TRUE.
            ENDIF
         ENDDO

         call read_par(kpp_3d_fields,ncid_landsea,'max_depth',1,1,
     +        kpp_3d_fields%ocdepth)

         WRITE(6,*) 'Read landsea mask'

         CALL close_landseafile
      ELSE
         DO ipt=1,npts
            kpp_3d_fields%L_OCEAN(ipt)=.FALSE.
         ENDDO
      ENDIF

      RETURN
      END
************************************************************
      SUBROUTINE init_cplwght
      
      IMPLICIT NONE
      INTEGER nuout,nuerr,start(2),count(2)
      INTEGER ix,jy,ipoint,cplwght_varid,status
      INTEGER ipoint_globe
      
      PARAMETER (nuout=6,nuerr=0)
#include <netcdf.inc>
#include <parameter.inc>
#include <location.com>
#include <couple.com>

      REAL*4 ixx, jyy, cplwght_in(NX_GLOBE,NY_GLOBE)
c      
c     If L_CPLWGHT has been set, then we will use the
c     NetCDF file to set values of cplwght over the
c     entire globe.
c
c     Otherwise, we will set the values ourselves, based
c     on the coupling region.
c
c     NPK 10/9/07 - R1
c
c      WRITE(6,*) 'L_CPLWGHT =',L_CPLWGHT
      IF (L_CPLWGHT) THEN
         CALL init_cplwghtfile
         start(1) = 1
         start(2) = 1
         count(1) = NX_GLOBE
         count(2) = NY_GLOBE
         WRITE(6,*) 'Reading coupling weight (alpha)'
         status=NF_INQ_VARID(ncid_cplwght,'alpha',cplwght_varid)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         status=NF_GET_VARA_REAL(ncid_cplwght,cplwght_varid,start,
     &        count,cplwght_in)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         WRITE(6,*) 'Coupling weight successfully read'
         DO ix=1,NX_GLOBE
            DO jy=1,NY_GLOBE
               ipoint_globe=(jy-1)*NX_GLOBE+ix
               cplwght(ipoint_globe)=cplwght_in(ix,jy)
            ENDDO
         ENDDO
         CALL close_cplwghtfile
      ELSE
         cplwght(:) = 2.   
      ENDIF
      DO ix=1,NX_GLOBE
         DO jy=1,NY_GLOBE
            ipoint_globe=(jy-1)*NX_GLOBE+ix            
            ixx=MIN(ix-ifirst,ilast-ix)
            jyy=MIN(jy-jfirst,jlast-jy)
            IF (ixx .GE. 0 .AND. jyy .GE. 0) THEN
c     Point is inside coupling domain.  
c     Set cplwght equal to one (if not already set from NetCDF file) 
c     to obtain model SSTs.
               cplwght(ipoint_globe) = MIN(cplwght(ipoint_globe),1.)
c               WRITE(6,*) 'ix=',ix,'jy=',jy,
c     +              'cplwght=',cplwght(ipoint_globe)
            ELSE
c     Point is outside coupling domain.
c     Set cplwght equal to a negative value to obtain
c     climatological SSTs (or persisted SSTs, IF (.NOT. L_UPDCLIM))
               cplwght(ipoint_globe) = -1.
            ENDIF
         ENDDO
      ENDDO
               
      RETURN
      
      END SUBROUTINE init_cplwght
