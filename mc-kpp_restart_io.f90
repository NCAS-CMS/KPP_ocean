SUBROUTINE MCKPP_RESTART_IO_READ(kpp_3d_fields,kpp_const_fields,restart_infile)
  
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
  ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
  
  ! Inputs 
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  CHARACTER(LEN=17) :: restart_infile
  
  WRITE(6,*) 'Total number of points = ',REAL(NPTS)*REAL(NZP1)
  IF ( REAL(NPTS)*REAL(NZP1) .LT. 3000000. ) THEN   
     OPEN(30,FILE=restart_infile,status='unknown',form='unformatted')
     READ(30) kpp_const_fields%time,kpp_3d_fields%U,kpp_3d_fields%X,kpp_3d_fields%CP,&
          kpp_3d_fields%rho,kpp_3d_fields%hmix,kpp_3d_fields%kmix,kpp_3d_fields%Sref,&
          kpp_3d_fields%SSref,kpp_3d_fields%Ssurf,kpp_3d_fields%Tref,kpp_3d_fields%old,&
          kpp_3d_fields%new,kpp_3d_fields%Us,kpp_3d_fields%Xs,kpp_3d_fields%hmixd
     CLOSE(30)
  ELSE
     OPEN(30,FILE=restart_infile//'.1',status='unknown',form='unformatted')
     READ(30) kpp_const_fields%time,kpp_3d_fields%U,kpp_3d_fields%X,kpp_3d_fields%CP,&
          kpp_3d_fields%rho,kpp_3d_fields%hmix,kpp_3d_fields%kmix,kpp_3d_fields%Sref,&
          kpp_3d_fields%SSref,kpp_3d_fields%Ssurf,kpp_3d_fields%Tref,kpp_3d_fields%old,&
          kpp_3d_fields%new
     CLOSE(30)
     OPEN(31,FILE=restart_infile//'.2',status='unknown',form='unformatted')
     READ(31) kpp_3d_fields%Us,kpp_3d_fields%Xs,kpp_3d_fields%hmixd
     CLOSE(31)
  ENDIF

  IF (abs(kpp_const_fields%time-kpp_const_fields%startt) .GT. 1.e-4) THEN 
     WRITE(nuerr,*) 'Start time doesn''t match the restart record'
     WRITE(nuerr,*) 'Start time in restart record = ',kpp_const_fields%time
     WRITE(nuerr,*) 'Start time in namelist = ',kpp_const_fields%startt
     CALL MCKPP_ABORT
  ENDIF
  
  RETURN
END SUBROUTINE MCKPP_RESTART_IO_READ

SUBROUTINE MCKPP_RESTART_IO_WRITE(kpp_3d_fields,kpp_const_fields,restart_outfile)
      
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
  ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>

  ! Inputs
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  CHARACTER(LEN=17) :: restart_outfile
  
  ! When the number of points in the model (NX*NY*NZP1) becomes
  ! quite large, we exceed the maximum size for Fortran unformatted
  ! binary files on certain machines.  The IF test below
  ! works around this by splitting the restart file in two.
  ! %Us and %Xs are the largest fields, so they get their own file.

  WRITE(6,*) 'Total number of points = ',REAL(NPTS)*REAL(NZP1)
  IF ( REAL(NPTS)*REAL(NZP1) .LT. 3000000. ) THEN
     OPEN(31,FILE=restart_outfile,status='unknown',form='unformatted')
     WRITE(31) kpp_const_fields%time,kpp_3d_fields%U,kpp_3d_fields%X,kpp_3d_fields%CP,&
          kpp_3d_fields%rho,kpp_3d_fields%hmix,kpp_3d_fields%kmix,kpp_3d_fields%Sref,&
          kpp_3d_fields%SSref,kpp_3d_fields%Ssurf,kpp_3d_fields%Tref,kpp_3d_fields%old,&
          kpp_3d_fields%new,kpp_3d_fields%Us,kpp_3d_fields%Xs,kpp_3d_fields%hmixd
     CLOSE(31)
  ELSE
     OPEN(31,FILE=restart_outfile//'.1',status='unknown',form='unformatted')
     WRITE(31) kpp_const_fields%time,kpp_3d_fields%U,kpp_3d_fields%X,kpp_3d_fields%CP,&
          kpp_3d_fields%rho,kpp_3d_fields%hmix,kpp_3d_fields%kmix,kpp_3d_fields%Sref,&
          kpp_3d_fields%SSref,kpp_3d_fields%Ssurf,kpp_3d_fields%Tref,kpp_3d_fields%old,&
          kpp_3d_fields%new
     CLOSE(31)
     OPEN(32,FILE=restart_outfile//'.2',status='unknown',form='unformatted')
     WRITE(32) kpp_3d_fields%Us,kpp_3d_fields%Xs,kpp_3d_fields%hmixd
     CLOSE(32)
  ENDIF
  
  RETURN
END SUBROUTINE MCKPP_RESTART_IO_WRITE
