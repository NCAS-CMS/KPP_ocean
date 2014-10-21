SUBROUTINE mckpp_initialize_ocean_profiles(kpp_3d_fields,kpp_const_fields)

  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)

#include <netcdf.inc>  
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
#include <initialcon.com>
#include <ocn_paras.com>
#include <ocn_state.com>
#include <constants.com>
#include <vert_pgrid.com>

! local
  INTEGER n,i,ipt
  INTEGER status,ncid
  REAL*4, allocatable :: var_in(:,:,:),z_in(:),x_in(:),y_in(:)
  
  INTEGER varid,dimid
  INTEGER nz_in,nx_in,ny_in
  
  INTEGER ix,iy,count(3),start(3)
  INTEGER kin,k
  REAL deltaz,deltavar,offset_sst
  
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  
  IF ( L_INITDATA ) THEN
     allocate(var_in(NX,NY,200))
     allocate(z_in(200))
     allocate(x_in(NX_GLOBE))
     allocate(y_in(NY_GLOBE))
     
     status=NF_OPEN(initdata_file,0,ncid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)     
     status=NF_INQ_DIMID(ncid,'longitude',dimid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_DIMLEN (ncid,dimid,nx_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'longitude',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,x_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     ix=1
     DO WHILE (abs(x_in(ix)-kpp_3d_fields%dlon(1)) .GT. 1.e-3)
        ix=ix+1
        IF (ix .GE. nx_in) THEN
           write(nuerr,*) 'Error reading initial conditions'
           write(nuerr,*) 'Can''t find longitude ',kpp_3d_fields%dlon(1),' in range ',x_in(1),x_in(nx_in)
           CALL MCKPP_ABORT
        ENDIF
     ENDDO
     start(1)=ix
     count(1)=nx

     status=NF_INQ_DIMID(ncid,'latitude',dimid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_DIMLEN (ncid,dimid,ny_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'latitude',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,y_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     iy=1
     DO WHILE (abs(y_in(iy)-kpp_3d_fields%dlat(1)) .GT. 1.e-3)
        iy=iy+1
        IF (iy .GE. ny_in) THEN
           write(nuerr,*) 'Error reading initial conditions'
           write(nuerr,*) 'Can''t find latitude ',kpp_3d_fields%dlat(1),' in range ',y_in(1),y_in(ny_in)
           CALL MCKPP_ABORT
        ENDIF
     ENDDO
     start(2)=iy
     count(2)=ny
     
     status=NF_INQ_DIMID(ncid,'zvel',dimid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_DIMLEN (ncid,dimid,nz_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     start(3)=1
     count(3)=nz_in
     
     status=NF_INQ_VARID(ncid,'zvel',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,z_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'u',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
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
                    kpp_3d_fields%U(ipt,k,1)=var_in(ix,iy,kin)+deltavar*(kpp_const_fields%zm(k)-z_in(kin))/deltaz
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ELSE
        write(nuerr,*) 'You have to interpolate'
     ENDIF
     write(6,*) 'read_init interpolated u'
     status=NF_INQ_VARID(ncid,'v',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
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
                    kpp_3d_fields%U(ipt,k,2)=var_in(ix,iy,kin)+deltavar*(kpp_const_fields%zm(k)-z_in(kin))/deltaz
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ELSE
        write(nuerr,*) 'You have to interpolate'
     ENDIF
     write(6,*) 'read_init interpolated v'
     ! Save initial currents in case they are needed to reinitalise
     ! dodgy profiles (see resetting routines in steves_3d_ocn.f)
     kpp_3d_fields%U_init=kpp_3d_fields%U
     
     status=NF_INQ_DIMID(ncid,'ztemp',dimid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_DIMLEN (ncid,dimid,nz_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     start(3)=1
     count(3)=nz_in
     
     status=NF_INQ_VARID(ncid,'ztemp',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,z_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'temp',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
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
                    kpp_3d_fields%X(ipt,k,1)=var_in(ix,iy,kin)+deltavar*(kpp_const_fields%zm(k)-z_in(kin))/deltaz
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ELSE
        write(nuerr,*) 'You have to interpolate'
     ENDIF
     
     ! KPP requires temperatures in CELSIUS.  If initial conditions
     ! are in Kelvin, subtract 273.15
     offset_sst = 0.
     DO ix=1,nx
        DO iy=1,ny
           IF (var_in(ix,iy,1) .gt. 200 .and. var_in(ix,iy,1) .lt. 400) offset_sst = kpp_const_fields%TK0
        END DO
     END DO
     kpp_3d_fields%X(:,:,1) = kpp_3d_fields%X(:,:,1) - offset_sst
     
     status=NF_INQ_DIMID(ncid,'zsal',dimid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_DIMLEN (ncid,dimid,nz_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     start(3)=1
     count(3)=nz_in
     
     status=NF_INQ_VARID(ncid,'zsal',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,z_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'sal',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
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
                    kpp_3d_fields%X(ipt,k,2)=var_in(ix,iy,kin)+deltavar*(kpp_const_fields%zm(k)-z_in(kin))/deltaz
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ELSE
        write(nuerr,*) 'You have to interpolate'
     ENDIF
  ELSE     
     write(nuerr,*) "No code for L_INITDATA=.FALSE."
     CALL MCKPP_ABORT
  ENDIF

  ! Calculate and remove reference salinity      

  DO ipt=1,npts
     kpp_3d_fields%Sref(ipt)=(kpp_3d_fields%X(ipt,1,2)+kpp_3d_fields%X(ipt,nzp1,2))/2.
     kpp_3d_fields%Ssref(ipt)=kpp_3d_fields%Sref(ipt)
     do i=1,nzp1
        kpp_3d_fields%X(ipt,i,2)=kpp_3d_fields%X(ipt,i,2)-kpp_3d_fields%Sref(ipt)
     enddo
     ! Initial surface temp
     kpp_3d_fields%Tref(ipt) = kpp_3d_fields%X(ipt,1,1)
     IF (L_SSref) THEN
        kpp_3d_fields%Ssurf(ipt)=kpp_3d_fields%SSref(ipt)
     ELSE
        kpp_3d_fields%Ssurf(ipt)=kpp_3d_fields%X(ipt,1,2)+kpp_3d_fields%Sref(ipt)
     ENDIF
  ENDDO
  
  RETURN
END SUBROUTINE mckpp_initialize_ocean_profiles

SUBROUTINE mckpp_initialize_ocean_model(kpp_3d_fields,kpp_const_fields)
  ! Initialize ocean model:
  ! Set coefficients for tridiagonal matrix solver.
  ! Compute hmix and diffusivity profiles for initial profile.
  ! Prepare for first time step.
  
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)

  ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
  
  ! Input
  TYPE(kpp_const_type) :: kpp_const_fields
  
  ! Output
  TYPE(kpp_3d_type) :: kpp_3d_fields
  
  ! Local
  TYPE(kpp_1d_type) :: kpp_1d_fields
  real dzb(NZ)              ! diff. between grid-levels below z(j)
  integer k,kmix0,n,l,ipt
  real hmix0,deltaz

  ! Compute factors for coefficients of tridiagonal matrix elements.
  ! tri(0     ,1,.........) : dt/h(1) factor for rhs flux
  ! tri(k=1:NZ,0,.........) : dt/h(k)/ {dzb(k-1)=z(k-1)-z(k)=dzabove}
  ! tri(k=1:NZ,1,.........) : dt/h(k)/ {dzb(k  )=z(k)-z(k+1)=dzbelow}

  DO k=1,NZ
     dzb(k) = kpp_const_fields%zm(k) - kpp_const_fields%zm(k+1)
  ENDDO
  
  kpp_const_fields%tri(0,1,1) = kpp_const_fields%dto/kpp_const_fields%hm(1)
  kpp_const_fields%tri(1,1,1) = kpp_const_fields%dto/kpp_const_fields%hm(1)/dzb(1)
  DO k=2,NZ
     kpp_const_fields%tri(k,1,1) = kpp_const_fields%dto/kpp_const_fields%hm(k)/dzb(k)
     kpp_const_fields%tri(k,0,1) = kpp_const_fields%dto/kpp_const_fields%hm(k)/dzb(k-1)
  ENDDO
  
  IF ( .NOT. kpp_const_fields%L_RESTART) THEN
    
     ! Determine hmix for initial profile:
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(kpp_3d_fields,kpp_const_fields)
!$OMP DO SCHEDULE(dynamic)
#endif
     DO ipt=1,npts
        CALL mckpp_fields_3dto1d(kpp_3d_fields,ipt,kpp_1d_fields)            
        IF (kpp_1d_fields%L_OCEAN) THEN
           kpp_1d_fields%L_INITFLAG=.TRUE.
           CALL MCKPP_PHYSICS_VERTICALMIXING(kpp_1d_fields,kpp_const_fields,hmix0,kmix0)
           kpp_1d_fields%L_INITFLAG=.FALSE.
           kpp_1d_fields%hmix = hmix0
           kpp_1d_fields%kmix = kmix0
           kpp_1d_fields%Tref = kpp_1d_fields%X(1,1)
           ! Evaluate initial fluxes (to write to output data file)
           DO k=1,NZ
              deltaz = 0.5*(kpp_const_fields%hm(k)+kpp_const_fields%hm(k+1))
              DO n=1,NSCLR
                 kpp_1d_fields%wX(k,n)=-kpp_1d_fields%difs(k)*&
                      ((kpp_1d_fields%X(k,n)-kpp_1d_fields%X(k+1,n))/deltaz-&
                      kpp_1d_fields%ghat(k)*kpp_1d_fields%wX(0,n))
              ENDDO
              IF(kpp_const_fields%LDD) kpp_1d_fields%wX(k,1)=-kpp_1d_fields%dift(k)*&
                   ((kpp_1d_fields%X(k,1)-kpp_1d_fields%X(k+1,1))/deltaz-kpp_1d_fields%ghat(k)*kpp_1d_fields%wX(0,1))
              kpp_1d_fields%wX(k,nsp1)= kpp_const_fields%grav * (kpp_1d_fields%talpha(k)*kpp_1d_fields%wX(k,1) - &
                   kpp_1d_fields%sbeta(k) * kpp_1d_fields%wX(k,2))
              DO  n=1,NVEL
                 kpp_1d_fields%wU(k,n)= -kpp_1d_fields%difm(k)*&
                      (kpp_1d_fields%U(k,n)-kpp_1d_fields%U(k+1,n))/deltaz
              ENDDO
           ENDDO
           
           ! Prepare for first time step
           
           ! indices for extrapolation
           kpp_1d_fields%old = 0
           kpp_1d_fields%new = 1               
           ! initialize array for extrapolating hmixd,Us,Xs
           kpp_1d_fields%hmixd(0) = kpp_1d_fields%hmix
           kpp_1d_fields%hmixd(1) = kpp_1d_fields%hmix
           DO k=1,NZP1
              DO l=1,NVEL
                 kpp_1d_fields%Us(k,l,0)=kpp_1d_fields%U(k,l)
                 kpp_1d_fields%Us(k,l,1)=kpp_1d_fields%U(k,l)
              ENDDO
              DO l=1,NSCLR
                 kpp_1d_fields%Xs(k,l,0)=kpp_1d_fields%X(k,l)
                 kpp_1d_fields%Xs(k,l,1)=kpp_1d_fields%X(k,l)
              ENDDO
           ENDDO
        ENDIF
        CALL mckpp_fields_1dto3d(kpp_1d_fields,ipt,kpp_3d_fields)
     ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif      
  ENDIF
  RETURN
END SUBROUTINE mckpp_initialize_ocean_model
