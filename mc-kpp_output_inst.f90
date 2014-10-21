SUBROUTINE mckpp_output_inst(kpp_3d_fields,kpp_const_fields,diag_num)      
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
! #include <output.com>
#include <times.com>
#include <ocn_advec.com>
#include <landsea.com>
#include <relax_3d.com>
#include <couple.com>
#include <fcorr_in.com>
#include <sfcorr_in.com>
  
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  REAL*4, allocatable :: varout(:,:,:), singout(:,:),temp_2d(:,:),temp_1d(:),temp_zprof(:,:)
  REAL*4, parameter :: missval=1.e20
  INTEGER start(4),count(4),status
  INTEGER i,j,k,ivar,diag_num
  INTEGER ix,iy,ipt,zprof
  
  count(1)=NX
  count(2)=NY
  start(1)=1
  start(2)=1
        
  write(nuout,*) 'Writing output for diagnostic ',diag_num,&
       'to start=',start,'count=',count,'zprof=',zprof
  
  IF (diag_num .le. N_VAROUTS) THEN 
     zprof=kpp_const_fields%zprof_varout_inst(diag_num)
     start(3)=1
     start(4)=kpp_const_fields%ntout_vec_inst(diag_num)
     count(3)=kpp_const_fields%zprofs_nvalid(zprof)
     count(4)=1
     allocate(varout(NX,NY,kpp_const_fields%zprofs_nvalid(zprof)))
     allocate(temp_2d(NPTS,NZP1))
     SELECT CASE (diag_num)
     CASE (1)
        temp_2d(:,:)=kpp_3d_fields%U(:,:,1)
     CASE (2)
        temp_2d(:,:)=kpp_3d_fields%U(:,:,2)
     CASE (3)
        temp_2d(:,:)=kpp_3d_fields%X(:,:,1)
     CASE (4)
        DO k=1,NZP1
           temp_2d(:,k)=kpp_3d_fields%X(:,k,2)+kpp_3d_fields%Sref(:)
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
     END SELECT
     
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

     status=NF_PUT_VARA_REAL(kpp_const_fields%ncid_out,kpp_const_fields%varid_vec_inst(diag_num),&
          start,count,VAROUT)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     kpp_const_fields%ntout_vec_inst(diag_num)=kpp_const_fields%ntout_vec_inst(diag_num)+1
  ELSE
     allocate(singout(NX,NY))      
     allocate(temp_1d(NPTS))
     start(3)=kpp_const_fields%ntout_sing_inst(diag_num-N_VAROUTS)
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
              temp_1d((iy-jfirst)*NX+ix-ifirst+1)=kpp_3d_fields%cplwght(ipt)                     
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
        WRITE(6,*) 'You need to add more outputs in OUTPUT_INST'
     END SELECT
     
     CALL MCKPP_REFORMAT_MASK_OUTPUT_1D(temp_1d,kpp_3d_fields%L_OCEAN,missval,singout)
     
     status=NF_PUT_VARA_REAL(kpp_const_fields%ncid_out,kpp_const_fields%varid_sing_inst(diag_num-N_VAROUTS),&
          start,count,SINGOUT)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     kpp_const_fields%ntout_sing_inst(diag_num-N_VAROUTS)=&
          kpp_const_fields%ntout_sing_inst(diag_num-N_VAROUTS)+1
  ENDIF
  start(3)=1
  
  RETURN
END SUBROUTINE mckpp_output_inst
