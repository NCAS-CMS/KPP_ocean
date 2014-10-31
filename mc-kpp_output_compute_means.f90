SUBROUTINE mckpp_output_compute_means(kpp_3d_fields,kpp_const_fields)
  
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)

  ! Automatically includes parameter.inc
#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  REAL, allocatable :: field(:,:),vec(:)
  INTEGER :: i,j,k,ivar,ix,iy,ipt,ipt_globe
  
  allocate(field(NPTS,NZP1))
  allocate(vec(NPTS))

  i=1
  DO ivar=1,N_VAROUTS
     IF (kpp_const_fields%ndt_varout_mean(ivar) .gt. 0) THEN
        !            WRITE(6,*) 'Computing means for ivar = ',ivar,'i=',i
        !            WRITE(6,*) 'ndt_varout_mean(ivar)=',ndt_varout_mean(ivar)
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
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,kpp_const_fields) &
!$OMP SHARED(i,field,ivar)
!$OMP DO SCHEDULE(static)
#endif
        DO j=1,NPTS
           IF (kpp_3d_fields%L_OCEAN(j)) THEN
              DO k=1,NZP1
                 kpp_3d_fields%VEC_mean(j,k,i)=&
                      field(j,k)/kpp_const_fields%ndt_varout_mean(ivar) + kpp_3d_fields%VEC_mean(j,k,i)
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
           IF (kpp_3d_fields%L_OCEAN(j)) &
                kpp_3d_fields%SCLR_mean(j,i)=vec(j)/kpp_const_fields%ndt_singout_mean(ivar)+&
                kpp_3d_fields%SCLR_mean(j,i)
        ENDDO
        i=i+1
     ENDIF
  ENDDO
  deallocate(field)
  deallocate(vec)
  
  RETURN
END SUBROUTINE mckpp_output_compute_means
