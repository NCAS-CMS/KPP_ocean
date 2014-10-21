#ifdef COUPLE
#ifdef OASIS3
SUBROUTINE MCKPP_COUPLING_OASIS3_INPUT(kpp_3d_fields,kpp_const_fields)
 
  ! Facilitate the exchange of coupled fields via the OASIS3 coupler.
  ! NPK 18/09/09 for the OASIS3 toy model - completed 28/09/09 - R3
  
  ! Use OASIS3 (PRISM) modules - the compiler must be able to find
  ! these via the Makefile.
  USE mod_kinds_model
  USE mod_prism_proto
  USE mod_prism_def_partition_proto
  USE mod_prism_put_proto
  USE mod_prism_get_proto
  USE mod_prism_grids_writing
  
  IMPLICIT NONE
  
  ! Standard include files
#include <netcdf.inc>
  ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
#include <kpp_oasis3.inc>
  
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  ! Local variables
  
  ! Note: "time_in_seconds" must be defined as INTEGER to
  ! be consistent with the definition in OASIS3.

  ! Note: the "temporary" variable must be defined
  ! as the OASIS3 type "ip_realwp_p".  OASIS3 must be compiled
  ! with default double precision (-fdefault-real-8 or similar).
#ifdef TOYCLIM
  REAL(KIND=ip_realwp_p) temporary(NPTS_GLOBE)
#else
  REAL(KIND=ip_realwp_p) temporary(NX_GLOBE,NY_GLOBE)
#endif
  REAL rain(NPTS),evap(NPTS)
  INTEGER i,j,ierror
  INTEGER time_in_seconds
  
  ! OASIS3 expects the time_in_seconds
  ! to be the time since the start of *this particular run*,
  ! not the time since the beginning of the initial run.
  ! NPK 21/12/09
  time_in_seconds=kpp_const_fields%spd*(kpp_const_fields%time-kpp_const_fields%startt)
  
  ! Get the coupled fields from the OASIS coupler.  Note that you
  ! can discard any fields you do not want by simply not defining
  ! a CASE for them.

  DO i=1,jpfldin
     CALL prism_get_proto(il_var_id_in(i),time_in_seconds,temporary,ierror)
     IF (ierror.NE.PRISM_Ok .and. ierror .LT. PRISM_Recvd) THEN
        WRITE(il_mparout) 'KPP: Received error from PRISM_Get_Proto =',ierror,&
             ' receiving variable ',cl_read(i),' at model time ',time_in_seconds,' sec.'
        WRITE(il_mparout) 'KPP: Aborting coupled integration ...'
        CALL prism_abort_proto(il_comp_id,'couple_io_oasis3.f','get')
     ELSE
        SELECT CASE (cl_read(i))
           ! For each field that we are coupling, use <TWOD_GLOBAL_ONED_REGIONAL>
           ! to transform that field from the global atmospheric grid to the
           ! KPP regional grid.
        CASE('HEATFLUX')
#ifdef TOYCLIM
           CALL MCKPP_ONED_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,kpp_3d_fields%lwf)
#else
           CALL MCKPP_TWOD_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,kpp_3d_fields%lwf)
#endif
        CASE('SOLAR')
#ifdef TOYCLIM
           CALL MCKPP_ONED_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,kpp_3d_fields%swf)
#else
           CALL MCKPP_TWOD_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,kpp_3d_fields%swf)
#endif
        CASE('PEN_SOL')
#ifdef TOYCLIM
           CALL MCKPP_ONED_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,kpp_3d_fields%swf)
#else
           CALL MCKPP_TWOD_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,kpp_3d_fields%swf)
#endif
           
        CASE('TAUX')
#ifdef TOYCLIM
           CALL MCKPP_ONED_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,kpp_3d_fields%taux)
#else
           CALL MCKPP_TWOD_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,kpp_3d_fields%taux)
#endif
        CASE('TAUY')
#ifdef TOYCLIM
           CALL MCKPP_ONED_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,kpp_3d_fields%tauy)
#else
           CALL MCKPP_TWOD_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,kpp_3d_fields%tauy)
#endif
        CASE ('TRAIN')
#ifdef TOYCLIM
           CALL MCKPP_ONED_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,rain)
#else
           CALL MCKPP_TWOD_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,rain)
#endif
        CASE ('EVAP2D')
#ifdef TOYCLIM
           CALL MCKPP_ONED_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,evap)
#else
           CALL MCKPP_TWOD_GLOBAL_ONED_REGIONAL(kpp_const_fields,temporary,evap)
#endif
        CASE DEFAULT
           WRITE(il_mparout,*) 'KPP: Discarding field ',cl_read(i)
        END SELECT
     ENDIF
  ENDDO
  DO i=1,NPTS
     kpp_3d_fields%rain(i)=rain(i)-evap(i)
  ENDDO
  RETURN      
END SUBROUTINE MCKPP_COUPLING_OASIS3_INPUT
      
SUBROUTINE MCKPP_COUPLING_OASIS3_OUTPUT(kpp_3d_fields,kpp_const_fields)
  
  ! Facilitate the exchange of coupled fields via the OASIS3 coupler.
  ! NPK 18/09/09 for the OASIS3 toy model - completed 28/09/09 - R3

  ! Use OASIS3 (PRISM) modules - the compiler must be able to find
  ! these via the Makefile.
  USE mod_kinds_model
  USE mod_prism_proto
  USE mod_prism_def_partition_proto
  USE mod_prism_put_proto
  USE mod_prism_get_proto
  USE mod_prism_grids_writing
       
  IMPLICIT NONE
  
  ! Standard include files
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
#include <kpp_oasis3.inc>

  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)

  ! Pass in the regional SST from the calling routine
  ! (usually <MAIN>).

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  ! Global SST and ICE variables that will be exported to OASIS
  REAL SST(NPTS_GLOBE)
  REAL ICE(NPTS_GLOBE)
  REAL SNOWDEPTH(NPTS_GLOBE)
  REAL ICEDEPTH(NPTS_GLOBE)
  REAL SURF_CURR_X(NPTS_GLOBE)
  REAL SURF_CURR_Y(NPTS_GLOBE)
  
  ! Temporary variable to allow the use of a SELECT CASE block.
  ! Note that this must be defined as "ip_realwp_p" and OASIS3 must
  ! be compiled with default REAL*8 (-fdefault-real-8 or similar).
#ifdef TOYCLIM
  REAL(KIND=ip_realwp_p) temporary(NX_GLOBE*NY_GLOBE)
#else
  REAL(KIND=ip_realwp_p) temporary(NX_GLOBE,NY_GLOBE)
#endif

  ! Other local variables
  ! Note: "time_in_seconds" must be INTEGER to agree with the
  ! definition in OASIS3.
  INTEGER i,ix,jy,ipoint_globe,ipoint,ierror
  INTEGER time_in_seconds
  
  ! If this is a restart run, we must correct for the
  ! non-zero start time.  OASIS3 expects the time_in_seconds
  ! to be the time since the start of *this particular run*,
  ! not the time since the beginning of the initial run.   
  ! NPK 21/12/09
  time_in_seconds=kpp_const_fields%spd*(kpp_const_fields%time-kpp_const_fields%startt)         
  WRITE(nuout,*) 'KPP: Time for coupling out is ',time_in_seconds,' seconds'

  ! Use the coupling weight to modify the SSTs before passing back to the coupler.
  ! N.B. : We want to do this whether or not the user has specified
  ! coupling weights.  Cplwght is automatically set IF (.NOT. L_CPLWGHT)   
  WRITE(nuout,*) ' Entering loop to weight SSTs'
  DO ix=1,NX_GLOBE
     DO jy=1,NY_GLOBE
        ipoint_globe = (jy-1)*NX_GLOBE+ix
        IF (kpp_3d_fields%cplwght(ipoint_globe) .LT. -1e-10) THEN              
           ! Point is outside the coupling domain; set to SST climatology
           SST(ipoint_globe) = kpp_3d_fields%sst(ix,jy)
           IF (.NOT. kpp_const_fields%L_CLIMCURR) THEN
              SURF_CURR_X(ipoint_globe)=0.
              SURF_CURR_Y(ipoint_globe)=0.
           ELSE
              SURF_CURR_X(ipoint_globe)=kpp_3d_fields%usf(ix,jy)
              SURF_CURR_Y(ipoint_globe)=kpp_3d_fields%vsf(ix,jy)
           ENDIF
        ELSE
           ! Point is inside the coupling domain; set to weighted value
           ipoint=(jy-kpp_const_fields%jfirst)*nx+(ix-kpp_const_fields%ifirst)+1             
           SST(ipoint_globe) = kpp_3d_fields%X(ipoint,1,1)*&
                kpp_3d_fields%cplwght(ipoint_globe)+kpp_3d_fields%sst(ix,jy)*&
                (1-kpp_3d_fields%cplwght(ipoint_globe))
           IF (kpp_const_fields%L_COUPLE_CURRENTS .AND. .NOT. kpp_const_fields%L_CLIMCURR) THEN
              SURF_CURR_X(ipoint_globe)=kpp_3d_fields%U(ipoint,1,1)*&
                   kpp_3d_fields%cplwght(ipoint_globe)
              SURF_CURR_Y(ipoint_globe)=kpp_3d_fields%U(ipoint,1,2)*&
                   kpp_3d_fields%cplwght(ipoint_globe)
           ELSEIF (kpp_const_fields%L_COUPLE_CURRENTS .AND. kpp_const_fields%L_CLIMCURR) THEN
              SURF_CURR_X(ipoint_globe)=kpp_3d_fields%usf(ix,jy)
              SURF_CURR_Y(ipoint_globe)=kpp_3d_fields%vsf(ix,jy)
           ENDIF
        ENDIF
        ice(ipoint_globe)=kpp_3d_fields%iceconc(ix,jy)

        ! If the user does not provide a climatological ice depth,
        ! then set the ice depth to be 2 metres.  Note that we set
        ! the depth to be 2 * the ice concentration because HadGEM3-A
        ! divides the provided depth by the ice concentration to obtain
        ! the mean depth over the ice-covered portion of the gridbox, 
        ! assuming that the ocean model provides the mean depth over the
        ! the entire gridbox.
        !     
        ! The previous, erroneous behaviour of setting our icedepth
        ! variable to 2m can be restored by setting L_BAD_ICE_DEPTH.            
        ! NPK updated 25/6/14.
        IF (.NOT. kpp_const_fields%L_CLIM_ICE_DEPTH) THEN
           IF (ice(ipoint_globe) .GE. 0) THEN
              IF (.NOT. kpp_const_fields%L_BAD_ICE_DEPTH) THEN
                 icedepth(ipoint_globe)=2.00*ice(ipoint_globe)
              ELSE
                 icedepth(ipoint_globe)=2.
              ENDIF
           ELSE
              icedepth(ipoint_globe)=0.00
           ENDIF
        ELSE
           icedepth(ipoint_globe)=kpp_3d_fields%icedepth(ix,jy)
        ENDIF

!     If there are no data available for the amount of snow on the ice,
!     then assume that there isn't any.  This is, again, based on the
!     default behavior in HadGEM3-A when running with AMIP-2 sea ice.
!     NPK 16/12/09 - R3
        IF (.NOT. kpp_const_fields%L_CLIM_SNOW_ON_ICE) THEN
           snowdepth(ipoint_globe)=0.00
        ELSE
           snowdepth(ipoint_globe)=kpp_3d_fields%snowdepth(ix,jy)
        ENDIF
     ENDDO
  ENDDO
  IF (kpp_const_fields%L_OUTKELVIN) SST = SST+kpp_const_fields%TK0
  WRITE(il_mparout,*) 'KPP: Finished creating coupling outputs'
  WRITE(nuout,*) 'KPP: Finished creating coupling outputs'

  DO i=1,jpfldout
     
     ! Export each field to OASIS.  Use a SELECT CASE block to avoid
     ! repeated bits of code.  Use the "temporary" variable to transfer
     ! the SST and ICE fields to the OASIS "ip_realwp_p" TYPE.
     SELECT CASE (cl_writ(i))
     CASE('OCN_SST')
#ifdef TOYCLIM
        temporary=SST
#else
        CALL MCKPP_ONED_GLOBAL_TWOD_GLOBAL(SST,temporary)
#endif
     CASE('OFRZN01')
#ifdef TOYCLIM
        temporary=ICE
#else
        CALL MCKPP_ONED_GLOBAL_TWOD_GLOBAL(ICE,temporary)
#endif
     CASE ('OSNWTN01')
#ifdef TOYCLIM
        temporary=SNOWDEPTH
#else
        CALL MCKPP_ONED_GLOBAL_TWOD_GLOBAL(SNOWDEPTH,temporary)
#endif            
     CASE ('OHICN01')
#ifdef TOYCLIM
        temporary=ICEDEPTH
#else
        CALL MCKPP_ONED_GLOBAL_TWOD_GLOBAL(ICEDEPTH,temporary)
#endif
     CASE ('SUNOCEAN')
#ifdef TOYCLIM
        temporary=SURF_CURR_X
#else
        IF (kpp_const_fields%L_COUPLE_CURRENTS) THEN
           CALL MCKPP_ONED_GLOBAL_TWOD_GLOBAL(SURF_CURR_X,temporary)
        ELSE
           temporary=0.
        ENDIF
#endif
     CASE ('SVNOCEAN')
#ifdef TOYCLIM
        temporary=SURF_CURR_Y
#else
        IF (kpp_const_fields%L_COUPLE_CURRENTS) THEN
           CALL MCKPP_ONED_GLOBAL_TWOD_GLOBAL(SURF_CURR_Y,temporary)
        ELSE
           temporary=0.
        ENDIF
#endif
     CASE DEFAULT
        ! The user should never see this - if they do, it means that the model
        ! is exporting a field that is not defined here (and likely not defined
        ! in the "namcouple" file either).
        WRITE(nuout,*) 'Unexpected CASE DEFAULT for i=',i
     END SELECT
     WRITE(nuout,*) 'KPP: Calling PRISM_Put_Proto for variable ',cl_writ(i)         
     CALL prism_put_proto(il_var_id_out(i),time_in_seconds,temporary,ierror)
     IF (ierror.NE.PRISM_Ok.and.ierror.LT.PRISM_Sent) THEN
        WRITE(nuout,*) 'KPP: Received error from PRISM_Put_Proto =',ierror,' sending variable ',&
             cl_writ(i),' at model time ',time_in_seconds,' sec.'
        WRITE(nuout,*) 'KPP: Aborting coupled integration ...'
        CALL prism_abort_proto(il_comp_id,'couple_io_oasis3.f','send')
     ELSE
        WRITE(nuout,*) 'KPP: Successfully called PRISM_Put_Proto for variable ',cl_writ(i)
     ENDIF
  ENDDO
 
  RETURN
END SUBROUTINE MCKPP_COUPLING_OASIS3_OUTPUT

SUBROUTINE MCKPP_COUPLING_OASIS3_FINALIZE
     
  ! Allows KPP to terminate its role in the coupled integration
  ! Essentially a wrapper around <prism_terminate_proto>
  ! NPK 28/09/09 - R3
     
  ! Use OASIS3 (PRISM) modules - the compiler must be able to find
  ! these via the Makefile.
  
  USE mod_kinds_model 
  USE mod_prism_proto
  
  IMPLICIT NONE

#include <kpp_oasis3.inc>

  ! Local variables
  INTEGER ierror
  
  ! Terminate the integration
  WRITE(6,*) 'KPP : Calling prism_terminate_proto(ierror)'
  CALL prism_terminate_proto(ierror)
  WRITE(6,*) 'KPP : Called prism_terminate_proto(ierror)'
  IF (ierror .NE. PRISM_Ok) THEN
     WRITE(il_mparout,*) 'KPP: Received error from PRISM_Terminate_Proto =',ierror,' when terminating model.'
  ENDIF
  WRITE(il_mparout,*) '---- End of the KPP integration ----'
  CLOSE(il_mparout)
    
  RETURN
END SUBROUTINE MCKPP_COUPLING_OASIS3_FINALIZE

SUBROUTINE MCKPP_ONED_GLOBAL_ONED_REGIONAL(kpp_const_fields,global,regional)

! Transforms a one-dimensional global field to a one-dimensional
! regional field.
! NPK 19/9/09 - R3 
  
  IMPLICIT NONE
#include <mc-kpp_3d_type.com>
!#include <couple.com>
  
  TYPE(kpp_const_type) :: kpp_const_fields
  REAL global(NPTS_GLOBE),regional(NPTS)
  INTEGER*4 ix,jy,ipoint,ipoint_globe
  
  DO jy=kpp_const_fields%jfirst,kpp_const_fields%jlast
     DO ix=kpp_const_fields%ifirst,kpp_const_fields%ilast
        ipoint_globe=jy*NX_GLOBE+ix
        ipoint=(jy-kpp_const_fields%jfirst)*nx+(ix-kpp_const_fields%ifirst)+1
        regional(ipoint)=global(ipoint_globe)
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE MCKPP_ONED_GLOBAL_ONED_REGIONAL

SUBROUTINE MCKPP_TWOD_GLOBAL_ONED_REGIONAL(kpp_const_fields,global,regional)

  ! Transforms a two-dimensional global field (e.g., one received
  ! from HadGEM3-A) to a one-dimensional regional field (as required for KPP).
  ! NPK 28/9/09 - R3

  IMPLICIT NONE
#include <mc-kpp_3d_type.com>
      
  TYPE(kpp_const_type) :: kpp_const_fields
  REAL global(NX_GLOBE,NY_GLOBE),regional(NPTS)
  INTEGER*4 ix,jy,ipoint,ipoint_globe

  DO jy=kpp_const_fields%jfirst,kpp_const_fields%jlast
     DO ix=kpp_const_fields%ifirst,kpp_const_fields%ilast
        ipoint=(jy-kpp_const_fields%jfirst)*NX+(ix-kpp_const_fields%ifirst)+1
        regional(ipoint)=global(ix,jy)
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE MCKPP_TWOD_GLOBAL_ONED_REGIONAL

SUBROUTINE MCKPP_ONED_GLOBAL_TWOD_GLOBAL(oned,twod)
  IMPLICIT NONE
      
! Transforms a one-dimensional global field (as produced by
! the combination of a KPP regional field and a global climatology)
! to a two-dimensional global field (as required by HadGEM3-A).
! NPK 28/9/09 - R3

#include <parameter.inc>
  REAL oned(NX_GLOBE*NY_GLOBE),twod(NX_GLOBE,NY_GLOBE)
  INTEGER*4 ix,jy,ipoint_globe
  
  DO jy=1,NY_GLOBE
     DO ix=1,NX_GLOBE
        ipoint_globe=(jy-1)*NX_GLOBE+ix
        twod(ix,jy)=oned(ipoint_globe)
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE MCKPP_ONED_GLOBAL_TWOD_GLOBAL

#endif /*OASIS3*/
#endif /*COUPLE*/
