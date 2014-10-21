SUBROUTINE mckpp_initialize_couplingweight(kpp_3d_fields)
  
  IMPLICIT NONE
  INTEGER nuout,nuerr,start(2),count(2)
  INTEGER ix,jy,ipoint,cplwght_varid,status
  INTEGER ipoint_globe
  
  PARAMETER (nuout=6,nuerr=0)
#include <mc-kpp_3d_type.com>
#include <netcdf.inc>
#include <couple.com>
  
  TYPE(kpp_3d_type) :: kpp_3d_fields
  REAL*4 ixx, jyy, cplwght_in(NX_GLOBE,NY_GLOBE)
  
!  If L_CPLWGHT has been set, then we will use the
!  NetCDF file to set values of cplwght over the
!  entire globe.
!  Otherwise, we will set the values ourselves, based
!  on the coupling region.
!  NPK 10/9/07 - R1
  
  IF (L_CPLWGHT) THEN
     status=NF_OPEN(cplwght_file,0,ncid_cplwght)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     start(1) = 1
     start(2) = 1
     count(1) = NX_GLOBE
     count(2) = NY_GLOBE
     WRITE(6,*) 'Reading coupling weight (alpha)'
     status=NF_INQ_VARID(ncid_cplwght,'alpha',cplwght_varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VARA_REAL(ncid_cplwght,cplwght_varid,start,count,cplwght_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     WRITE(6,*) 'Coupling weight successfully read'
     DO ix=1,NX_GLOBE
        DO jy=1,NY_GLOBE
           ipoint_globe=(jy-1)*NX_GLOBE+ix
           kpp_3d_fields%cplwght(ipoint_globe)=cplwght_in(ix,jy)
        ENDDO
     ENDDO
     status=NF_CLOSE(ncid_cplwght)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ELSE
     kpp_3d_fields%cplwght(:) = 2.   
  ENDIF
  DO ix=1,NX_GLOBE
     DO jy=1,NY_GLOBE
        ipoint_globe=(jy-1)*NX_GLOBE+ix            
        ixx=MIN(ix-ifirst,ilast-ix)
        jyy=MIN(jy-jfirst,jlast-jy)
        IF (ixx .GE. 0 .AND. jyy .GE. 0) THEN
           ! Point is inside coupling domain.  
           ! Set cplwght equal to one (if not already set from NetCDF file) 
           ! to obtain model SSTs.
           kpp_3d_fields%cplwght(ipoint_globe) = MIN(kpp_3d_fields%cplwght(ipoint_globe),1.)
           ipoint=NX*jyy+ixx
           IF (kpp_3d_fields%L_OCEAN(ipoint) .and. kpp_3d_fields%ocdepth(ipoint) .gt. 100) THEN
              kpp_3d_fields%L_OCEAN(ipoint)=.FALSE.
              kpp_3d_fields%cplwght(ipoint_globe)=0
              WRITE(6,*) 'Overwriting coupling mask at'
              WRITE(6,*) 'ixx=',ixx,'jyy=',jyy,'ipoint_globe=',ipoint_globe,'ipoint=',ipoint,'cplwght=',&
                   kpp_3d_fields%cplwght(ipoint_globe),'ocdepth=',kpp_3d_fields%ocdepth(ipoint),&
                   'L_OCEAN=',kpp_3d_fields%L_OCEAN(ipoint)
           ENDIF
        ELSE
           ! Point is outside coupling domain.
           ! Set cplwght equal to a negative value to obtain
           ! climatological SSTs (or persisted SSTs, IF (.NOT. L_UPDCLIM))
           kpp_3d_fields%cplwght(ipoint_globe) = -1.
        ENDIF
     ENDDO
  ENDDO
  
  RETURN  
END SUBROUTINE mckpp_initialize_couplingweight
