SUBROUTINE mckpp_output_compute_ranges(kpp_3d_fields,kpp_const_fields)

  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  ! Automatically includes parameter.inc
#include <mc-kpp_3d_type.com>
  
  REAL, allocatable :: field(:,:),vec(:)
  
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER i,j,k,ivar,ix,iy,ipt,ipt_globe
  
  allocate(field(NPTS,NZP1))
  allocate(vec(NPTS))

  i=1
  DO ivar=1,N_VAROUTS
     IF (kpp_const_fields%ndt_varout_range(ivar) .gt. 0) THEN
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
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,kpp_const_fields)&
!$OMP SHARED(i,field)
!$OMP DO SCHEDULE(static)
#endif
        DO j=1,NPTS
           IF (kpp_3d_fields%L_OCEAN(j)) THEN
              DO k=1,NZP1
                 IF (field(j,k) .lt. kpp_3d_fields%VEC_range(j,k,i,1)) &
                      kpp_3d_fields%VEC_range(j,k,i,1)=field(j,k)
                 IF (field(j,k) .gt. kpp_3d_fields%VEC_range(j,k,i,2)) &
                      kpp_3d_fields%VEC_range(j,k,i,2)=field(j,k)
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
     IF (kpp_const_fields%ndt_singout_mean(ivar) .gt. 0) THEN
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
           DO ix=kpp_const_fields%ifirst,kpp_const_fields%ilast
              DO iy=kpp_const_fields%jfirst,kpp_const_fields%jlast
                 ipt_globe=(iy-1)*NX_GLOBE+ix
                 ipt=(iy-kpp_const_fields%jfirst)*nx+(ix-kpp_const_fields%ifirst+1)
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
              IF (vec(j) .lt. kpp_3d_fields%SCLR_range(j,i,1)) &
                   kpp_3d_fields%SCLR_range(j,i,1)=vec(j)
              IF (vec(j) .gt. kpp_3d_fields%SCLR_range(j,i,2)) &
                   kpp_3d_fields%SCLR_range(j,i,2)=vec(j)
           ENDIF
        ENDDO
     ENDIF
     i=i+1      
  ENDDO
  deallocate(field)
  deallocate(vec)
  
  RETURN
END SUBROUTINE mckpp_output_compute_ranges
