SUBROUTINE mckpp_initialize_geography(lstretchgrid,dscale,kpp_3d_fields,kpp_const_fields)
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
  ! Automatically includes parameter.inc
#include <mc-kpp_3d_type.com>
#include <netcdf.inc>
#include <constants.com>
#include <vert_pgrid.com>

  ! Inputs
  LOGICAL lstretchgrid
  REAL dscale
  TYPE(kpp_3d_type)    :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  ! Local Variables
  REAL sumh,hsum,dfac,sk
  REAL*4 vgrid_in(NZ)
  INTEGER i,ipt,ncid,status,dimid,varid

  ! define vertical grid fields
  IF (L_VGRID_FILE) THEN 
     WRITE(6,*) 'Reading vertical grid from file ',vgrid_file
     status=NF_OPEN(vgrid_file,0,ncid)
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'d',varid)
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,vgrid_in)
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     kpp_const_fields%dm(1:NZ)=vgrid_in
     status=NF_OPEN(vgrid_file,0,ncid)
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'h',varid)
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,vgrid_in)
     kpp_const_fields%hm(1:NZ)=vgrid_in
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_OPEN(vgrid_file,0,ncid)
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'z',varid)
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,vgrid_in)
     kpp_const_fields%zm(1:NZ)=vgrid_in
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_CLOSE(ncid)
     DMAX=-1.*(kpp_const_fields%zm(NZ)-kpp_const_fields%hm(NZ))
  ELSE      
     IF (lstretchgrid) THEN
        sumh = 0.0
        dfac = 1.0 - exp(-dscale)
        DO i = 1,NZ
           sk = - (float(i)-0.5)/float(NZ)
           kpp_const_fields%hm(i) = DMAX*dfac/float(NZ)/dscale / ( 1.0 + sk*dfac )
           sumh = sumh + kpp_const_fields%hm(i)
        ENDDO
     ENDIF
        
     ! layer thickness h, layer grids zgrid, interface depths d
     hsum = 0.0
     kpp_const_fields%dm(0) = 0.0
     DO i=1,NZ
        if(lstretchgrid) then
           kpp_const_fields%hm(i) = kpp_const_fields%hm(i) * DMAX / sumh 
        else   
           kpp_const_fields%hm(i) = DMAX / real(NZ) 
        endif
        kpp_const_fields%zm(i) =  0.0 - (hsum + 0.5 * kpp_const_fields%hm(i) )
        hsum = hsum + kpp_const_fields%hm(i)
        kpp_const_fields%dm(i) = hsum
     ENDDO
  ENDIF
  kpp_const_fields%hm(nzp1) = 1.e-10 
  kpp_const_fields%zm(nzp1) = -DMAX
  
  ! Enforce minimum value of Coriolis parameter equal to 2.5 degrees latitude
  DO ipt=1,npts
     if(abs(kpp_3d_fields%dlat(ipt)).lt.2.5) then    
        kpp_3d_fields%f(ipt) = 2. * (twopi/86164.) * &
             sin(2.5*twopi/360.)*SIGN(1.,kpp_3d_fields%dlat(ipt))
     else  
        kpp_3d_fields%f(ipt) = 2. * (twopi/86164.) * &
             sin(kpp_3d_fields%dlat(ipt)*twopi/360.)
     endif
  ENDDO  
  
  RETURN
END SUBROUTINE mckpp_initialize_geography
