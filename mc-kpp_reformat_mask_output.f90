SUBROUTINE mckpp_reformat_mask_output_1d(oned_in,mask,missval,twod_out)
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
END SUBROUTINE mckpp_reformat_mask_output_1d

SUBROUTINE mckpp_reformat_mask_output_2d(twod_in,nz_in,mask,missval,threed_out)
  
  IMPLICIT NONE
#include <parameter.inc>
  
  INTEGER,intent(in) :: nz_in
  REAL*4,intent(in) :: twod_in(NPTS,nz_in),missval
  REAL*4,intent(out) :: threed_out(NX,NY,nz_in)
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
END SUBROUTINE mckpp_reformat_mask_output_2d
