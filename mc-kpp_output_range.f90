SUBROUTINE mckpp_output_range(kpp_3d_fields,kpp_const_fields,diag_num)
  IMPLICIT NONE

  INTEGER nuout,nuerr
  PARAMETER(nuout=6,nuerr=0)
  
#include <netcdf.inc>      
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
! #include 'output.com'
#include <landsea.com>
#include <times.com>
#include <timocn.com>
#include <ocn_advec.com>
#include <couple.com>
#include <constants.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  REAL*4,allocatable :: VAROUT(:,:,:), SINGOUT(:,:),temp_2d(:,:),temp_zprof(:,:)
  REAL*4, parameter :: missval=1.e20
  
  INTEGER i,j,ivar,ipt,ix,iy,start(4),count(4),k,status,diag_num,range_num,zprof
      
  count(1)=NX
  count(2)=NY
  start(1)=1
  start(2)=1
        
  allocate(SINGOUT(NX,NY))
  allocate(VAROUT(NX,NY,NZP1))
  allocate(temp_2d(NPTS,NZP1))
  DO j=1,2
     IF (diag_num .le. N_VAROUTS) THEN 
        zprof=kpp_const_fields%zprof_varout_range(diag_num)
        start(3)=1
        start(4)=kpp_const_fields%ntout_vec_range(diag_num)
        count(3)=kpp_const_fields%zprofs_nvalid(zprof)
        count(4)=1
        range_num=1
        DO i=1,diag_num-1
           IF (kpp_const_fields%ndt_varout_range(i) .gt. 0) range_num=range_num+1
        ENDDO
        SELECT CASE (diag_num)
        CASE (4) 
           DO k=1,NZP1
              temp_2d(:,k)=kpp_3d_fields%VEC_range(:,k,range_num,j)+kpp_3d_fields%Sref(:)
           ENDDO
        CASE (12,13,14)
           temp_2d(:,1)=0.
           temp_2d(:,2:NZP1)=kpp_3d_fields%VEC_range(:,2:NZP1,range_num,j)
        CASE (18,19)
           temp_2d(:,1:NZ)=kpp_3d_fields%VEC_range(:,1:NZ,range_num,j)
           temp_2d(:,NZP1)=0.               
        CASE DEFAULT
           temp_2d(:,:)=kpp_3d_fields%VEC_range(:,:,range_num,j)
        END SELECT
        
        IF (zprof .gt. 0) THEN
           k=1
           allocate(temp_zprof(NPTS,kpp_const_fields%zprofs_nvalid(zprof)))
           DO i=1,NZP1
              IF (kpp_const_fields%zprofs_mask(i,zprof)) THEN
                 temp_zprof(:,k)=temp_2d(:,k)
                 k=k+1
              ENDIF
           ENDDO
           CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_zprof,kpp_const_fields%zprofs_nvalid(zprof),&
                kpp_3d_fields%L_OCEAN,missval,varout)
           deallocate(temp_zprof)
        ELSE       
           CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_2d,NZP1,kpp_3d_fields%L_OCEAN,missval,varout)
        ENDIF
        SELECT CASE (j)
        CASE (1)                              
           WRITE(6,*) 'Writing to ncid=',kpp_const_fields%min_ncid_out,&
                ' varid=',kpp_const_fields%varid_vec_range(diag_num),&
                ' with start =',start,' and count =',count
           status=NF_PUT_VARA_REAL(kpp_const_fields%min_ncid_out,kpp_const_fields%varid_vec_range(diag_num),&
                start,count,VAROUT)
           IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
           kpp_3d_fields%VEC_range(:,:,range_num,j)=2e20
        CASE (2)           
           WRITE(6,*) 'Writing to ncid=',kpp_const_fields%max_ncid_out,&
                ' varid=',kpp_const_fields%varid_vec_range(diag_num),&
                ' with start =',start,' and count =',count
           status=NF_PUT_VARA_REAL(kpp_const_fields%max_ncid_out,kpp_const_fields%varid_vec_range(diag_num),&
                start,count,VAROUT)
           IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
           kpp_3d_fields%VEC_range(:,:,range_num,j)=2e20
        END SELECT        
        kpp_const_fields%ntout_vec_range(diag_num)=kpp_const_fields%ntout_vec_range(diag_num)+1
     ELSE
        start(3)=kpp_const_fields%ntout_sing_range(diag_num-N_VAROUTS)
        count(3)=1
        range_num=1
        DO i=1,(diag_num-N_VAROUTS-1)
           IF (kpp_const_fields%ndt_singout_range(i) .gt. 0) range_num=range_num+1
        ENDDO
        SINGOUT(:,:) = missval
        DO ix=1,nx
           DO iy=1,ny
              ipt=(iy-1)*nx+ix
              IF (kpp_3d_fields%L_OCEAN(ipt)) SINGOUT(ix,iy)=kpp_3d_fields%SCLR_range(ipt,range_num,j)
           ENDDO
        ENDDO
        SELECT CASE (j)
        CASE (1)               
           status=NF_PUT_VARA_REAL(kpp_const_fields%min_ncid_out,&
                kpp_const_fields%varid_sing_range(diag_num-N_VAROUTS),start,count,SINGOUT)
           IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
           kpp_3d_fields%SCLR_range(:,range_num,j)=2e20
        CASE (2)
           status=NF_PUT_VARA_REAL(kpp_const_fields%max_ncid_out,&
                kpp_const_fields%varid_sing_range(diag_num-N_VAROUTS),start,count,SINGOUT)
           IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
           kpp_3d_fields%SCLR_range(:,range_num,j)=-2e20
        END SELECT
        kpp_const_fields%ntout_sing_range(diag_num-N_VAROUTS)=&
             kpp_const_fields%ntout_sing_range(diag_num-N_VAROUTS)+1
     ENDIF
  ENDDO
  start(3)=1
  
  RETURN
END SUBROUTINE mckpp_output_range
