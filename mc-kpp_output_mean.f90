SUBROUTINE mckpp_output_mean(kpp_3d_fields,kpp_const_fields,diag_num)
  IMPLICIT NONE
      
  INTEGER nuout,nuerr
  PARAMETER(nuout=6,nuerr=0)

#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
! #include <output.com>
#include <landsea.com>
#include <times.com>
#include <timocn.com>
#include <ocn_advec.com>
#include <couple.com>
#include <constants.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  REAL*4,allocatable :: VAROUT(:,:,:), SINGOUT(:,:), temp_2d(:,:),temp_zprof(:,:)
  REAL*4, parameter :: missval=1.e20

  INTEGER i,j,ivar,ipt,ix,iy,start(4),count(4),k,status,diag_num,mean_num,zprof
  
  count(1)=NX
  count(2)=NY  
  start(1)=1
  start(2)=1
  
  IF (diag_num .le. N_VAROUTS) THEN
     zprof=kpp_const_fields%zprof_varout_mean(diag_num)
     start(3)=1
     start(4)=kpp_const_fields%ntout_vec_mean(diag_num)
     count(3)=kpp_const_fields%zprofs_nvalid(zprof)
     count(4)=1
     allocate(VAROUT(NX,NY,kpp_const_fields%zprofs_nvalid(zprof)))
     allocate(temp_2d(NPTS,NZP1))
     mean_num=1
     DO i=1,diag_num-1
        IF (kpp_const_fields%ndt_varout_mean(i) .gt. 0) mean_num=mean_num+1
     ENDDO
     SELECT CASE (diag_num)
     CASE (4) 
        DO k=1,NZP1
           temp_2d(:,k)=kpp_3d_fields%VEC_mean(:,k,mean_num)+kpp_3d_fields%Sref(:)
        ENDDO
     CASE (12,13,14)
        temp_2d(:,1)=0.
        temp_2d(:,2:NZP1)=kpp_3d_fields%VEC_mean(:,2:NZP1,mean_num)
     CASE (18,19)
        temp_2d(:,1:NZ)=kpp_3d_fields%VEC_mean(:,1:NZ,mean_num)
        temp_2d(:,NZP1)=0.               
     CASE DEFAULT
        temp_2d(:,:)=kpp_3d_fields%VEC_mean(:,:,mean_num)
     END SELECT
     WRITE(6,*) 'In WRITE_MEANS for diag_num=',diag_num,'zprof=',zprof
     
     IF (zprof .gt. 0) THEN
        j=1
        allocate(temp_zprof(NPTS,kpp_const_fields%zprofs_nvalid(zprof)))
        DO i=1,NZP1
           IF (kpp_const_fields%zprofs_mask(i,zprof)) THEN
              temp_zprof(:,j)=temp_2d(:,j)
              j=j+1
           ENDIF
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_zprof,kpp_const_fields%zprofs_nvalid(zprof),&
             kpp_3d_fields%L_OCEAN,missval,varout)
     ELSE       
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_2d,NZP1,kpp_3d_fields%L_OCEAN,missval,varout)
     ENDIF
     
     status=NF_PUT_VARA_REAL(kpp_const_fields%mean_ncid_out,kpp_const_fields%varid_vec_mean(diag_num),&
          start,count,VAROUT)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     kpp_3d_fields%VEC_mean(:,:,mean_num)=0.
     kpp_const_fields%ntout_vec_mean(diag_num)=kpp_const_fields%ntout_vec_mean(diag_num)+1
     i=i+1
  ELSE
     allocate(SINGOUT(NX,NY))     
     start(3)=kpp_const_fields%ntout_sing_mean(diag_num-N_VAROUTS)
     count(3)=1 
     mean_num=1
     DO i=1,(diag_num-N_VAROUTS-1)
        IF (kpp_const_fields%ndt_singout_mean(i) .gt. 0) mean_num=mean_num+1
     ENDDO
     SINGOUT(:,:) = missval
     DO ix=1,nx
        DO iy=1,ny
           ipt=(iy-1)*nx+ix
           IF (kpp_3d_fields%L_OCEAN(ipt)) SINGOUT(ix,iy)=kpp_3d_fields%SCLR_mean(ipt,mean_num)
        ENDDO
     ENDDO
     WRITE(nuout,*) 'In write_means for singout, i=',mean_num,'diag_num=',diag_num
     status=NF_PUT_VARA_REAL(kpp_const_fields%mean_ncid_out,&
          kpp_const_fields%varid_sing_mean(diag_num-N_VAROUTS),start,count,SINGOUT)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     kpp_const_fields%ntout_sing_mean(diag_num-N_VAROUTS)=&
          kpp_const_fields%ntout_sing_mean(diag_num-N_VAROUTS)+1
     kpp_3d_fields%SCLR_mean(:,mean_num)=0.
     i=i+1         
  ENDIF
  start(3)=1
  
  RETURN
END SUBROUTINE mckpp_output_mean
