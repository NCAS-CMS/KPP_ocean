#ifdef COUPLE
#ifdef OASIS3
      SUBROUTINE mpi1_oasis3_input(solar,non_solar,PminusE,ustress,
     +     vstress,kpp_const_fields)
c
c     Facilitate the exchange of coupled fields via the OASIS3 coupler.
c     NPK 18/09/09 for the OASIS3 toy model - completed 28/09/09 - R3
c
c     Use OASIS3 (PRISM) modules - the compiler must be able to find
c     these via the Makefile.
c
      USE mod_kinds_model
      USE mod_prism_proto
      USE mod_prism_def_partition_proto
      USE mod_prism_put_proto
      USE mod_prism_get_proto
      USE mod_prism_grids_writing
c
      IMPLICIT NONE
c
c     Standard include files
c
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
#include <kpp_oasis3.inc>
#include <times.com>
#include <constants.com>
      TYPE(kpp_const_type) :: kpp_const_fields
c
c     Output variables on the KPP regional grid - returned to
c     the calling routine (usually <fluxes>).
c
      REAL solar(NPTS)
      REAL non_solar(NPTS)
      REAL PminusE(NPTS)
      REAL ustress(NPTS)
      REAL vstress(NPTS)
c
c     Local variables
c     
c     Note: "time_in_seconds" must be defined as INTEGER to
c     be consistent with the definition in OASIS3.
c
c     Note: the "temporary" variable must be defined
c     as the OASIS3 type "ip_realwp_p".  OASIS3 must be compiled
c     with default double precision (-fdefault-real-8 or similar).
c
#ifdef TOYCLIM
      REAL(KIND=ip_realwp_p) temporary(NPTS_GLOBE)
#else
      REAL(KIND=ip_realwp_p) temporary(NX_GLOBE,NY_GLOBE)
#endif
      REAL rain(NPTS),evap(NPTS)
      INTEGER i,j,ierror
      INTEGER time_in_seconds
c
c     OASIS3 expects the time_in_seconds
c     to be the time since the start of *this particular run*,
c     not the time since the beginning of the initial run.
c     
c     NPK 21/12/09
c
      time_in_seconds=kpp_const_fields%spd*(kpp_const_fields%time-
     +     kpp_const_fields%startt)
                  
c      
c     Get the coupled fields from the OASIS coupler.  Note that you
c     can discard any fields you do not want by simply not defining
c     a CASE for them.
c
      DO i=1,jpfldin
         CALL prism_get_proto(il_var_id_in(i),
     +        time_in_seconds,temporary,ierror)
         IF (ierror.NE.PRISM_Ok .and. ierror .LT. PRISM_Recvd) THEN
            WRITE(il_mparout) 'KPP: Received error from ',
     +           'PRISM_Get_Proto =',ierror,' receiving variable ',
     +           cl_read(i),' at model time ',time_in_seconds,' sec.'
            WRITE(il_mparout) 'KPP: Aborting coupled integration ...'
            CALL prism_abort_proto(il_comp_id,'couple_io_oasis3.f',
     +           'get')
         ELSE
            SELECT CASE (cl_read(i))
c
c     For each field that we are coupling, use <TWOD_GLOBAL_ONED_REGIONAL>
c     to transform that field from the global atmospheric grid to the
c     KPP regional grid.
c
            CASE('HEATFLUX')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,non_solar)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,non_solar)
#endif
            CASE('SOLAR')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,solar)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,solar)
#endif
            CASE('TAUX')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,ustress)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,ustress)
#endif
            CASE('TAUY')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,vstress)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,vstress)
#endif
            CASE ('TRAIN')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,rain)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,rain)
#endif
            CASE ('EVAP2D')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,evap)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,evap)
#endif
            CASE DEFAULT
               WRITE(il_mparout,*) 'KPP: Discarding field ',
     +              cl_read(i)
            END SELECT
         ENDIF
      ENDDO
      DO i=1,NPTS
         PminusE(i)=rain(i)-evap(i)
      ENDDO

c      IF (READ_FROM_NETCDF)
c     +     CALL mpi1_oasis3_read_netcdf(solar,non_solar,PminusE,
c     +     ustress,vstress)
c      IF (WRITE_TO_NETCDF) 
c     +     CALL mpi1_oasis3_write_netcdf(solar,non_solar,PminusE,
c     +     ustress,vstress)         
c
c     End
c     
      RETURN      
      END SUBROUTINE mpi1_oasis3_input
      
      SUBROUTINE mpi1_oasis3_output(kpp_3d_fields,kpp_const_fields)
c     
c     Facilitate the exchange of coupled fields via the OASIS3 coupler.
c     NPK 18/09/09 for the OASIS3 toy model - completed 28/09/09 - R3
c
c     Use OASIS3 (PRISM) modules - the compiler must be able to find
c     these via the Makefile.
c
      USE mod_kinds_model
      USE mod_prism_proto
      USE mod_prism_def_partition_proto
      USE mod_prism_put_proto
      USE mod_prism_get_proto
      USE mod_prism_grids_writing
c      
      IMPLICIT NONE
c
c     Standard include files
c
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include <kpp_3d_type.com>
#include <kpp_oasis3.inc>
#include <times.com>
#include <couple.com>
c#include <location.com>
#include <constants.com>
#include <currclim.com>
c#include <initialcon.com>

      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
c
c     Pass in the regional SST from the calling routine
c     (usually <MAIN>).
c
      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
c
c     "SST_in" and "ICE_in" are common-block variables from
c     the netCDF input routines.
c
      REAL SST_in(NX_GLOBE,NY_GLOBE,1)
      REAL ICE_in(NX_GLOBE,NY_GLOBE,1)
      REAL ICEDEPTH_in(NX_GLOBE,NY_GLOBE,1)
      REAL SNOWDEPTH_in(NX_GLOBE,NY_GLOBE,1)
      REAL usf_in(NX_GLOBE,NY_GLOBE)
      REAL vsf_in(NX_GLOBE,NY_GLOBE)
c
c     Global SST and ICE variables that will be exported to OASIS
c
      REAL SST(NPTS_GLOBE)
      REAL ICE(NPTS_GLOBE)
      REAL SNOWDEPTH(NPTS_GLOBE)
      REAL ICEDEPTH(NPTS_GLOBE)
      REAL SURF_CURR_X(NPTS_GLOBE)
      REAL SURF_CURR_Y(NPTS_GLOBE)
c
c     Temporary variable to allow the use of a SELECT CASE block.
c     Note that this must be defined as "ip_realwp_p" and OASIS3 must
c     be compiled with default REAL*8 (-fdefault-real-8 or similar).
c
#ifdef TOYCLIM
      REAL(KIND=ip_realwp_p) temporary(NX_GLOBE*NY_GLOBE)
#else
      REAL(KIND=ip_realwp_p) temporary(NX_GLOBE,NY_GLOBE)
#endif
c
c     Other local variables
c     Note: "time_in_seconds" must be INTEGER to agree with the
c     definition in OASIS3.
c
      INTEGER i,ix,jy,ipoint_globe,ipoint,ierror
      INTEGER time_in_seconds
c
c     COMMON block for SST_in and ICE_in
c
      COMMON /save_sstin/ SST_in,ICE_in,icedepth_in,snowdepth_in,
     +     usf_in,vsf_in
c
c     If this is a restart run, we must correct for the
c     non-zero start time.  OASIS3 expects the time_in_seconds
c     to be the time since the start of *this particular run*,
c     not the time since the beginning of the initial run.
c     
c     NPK 21/12/09
c
      time_in_seconds=kpp_const_fields%spd*(kpp_const_fields%time-
     +     kpp_const_fields%startt)         
      WRITE(nuout,*) 'KPP: Time for coupling out is ',
     +     time_in_seconds,' seconds'
c
c     Use the coupling weight to modify the SSTs before passing
c     back to the coupler.
c     N.B. : We want to do this whether or not the user has specified
c     coupling weights.  Cplwght is automatically set IF (.NOT. L_CPLWGHT)
c      
      WRITE(nuout,*) ' Entering loop to weight SSTs'
      DO ix=1,NX_GLOBE
         DO jy=1,NY_GLOBE
            ipoint_globe = (jy-1)*NX_GLOBE+ix
            IF (kpp_3d_fields%cplwght(ipoint_globe) .LT. -1e-10) THEN              
c     Point is outside the coupling domain; set to SST climatology
               SST(ipoint_globe) = SST_in(ix,jy,1)
               IF (.NOT. L_CLIMCURR) THEN
                  SURF_CURR_X(ipoint_globe)=0.
                  SURF_CURR_Y(ipoint_globe)=0.
               ELSE
                  SURF_CURR_X(ipoint_globe)=usf_in(ix,jy)
                  SURF_CURR_Y(ipoint_globe)=vsf_in(ix,jy)
               ENDIF
            ELSE
c     Point is inside the coupling domain; set to weighted value
               ipoint=(jy-jfirst)*nx+(ix-ifirst)+1             
               SST(ipoint_globe) = kpp_3d_fields%X(ipoint,1,1)*
     +              kpp_3d_fields%cplwght(ipoint_globe)+SST_in(ix,jy,1)*
     +              (1-kpp_3d_fields%cplwght(ipoint_globe))
               IF (L_COUPLE_CURRENTS .AND. .NOT. L_CLIMCURR) THEN
                  SURF_CURR_X(ipoint_globe)=kpp_3d_fields%U(ipoint,1,1)*
     +                 kpp_3d_fields%cplwght(ipoint_globe)
                  SURF_CURR_Y(ipoint_globe)=kpp_3d_fields%U(ipoint,1,2)*
     +                 kpp_3d_fields%cplwght(ipoint_globe)
               ELSEIF (L_COUPLE_CURRENTS .AND. L_CLIMCURR) THEN
                  SURF_CURR_X(ipoint_globe)=usf_in(ix,jy)
                  SURF_CURR_Y(ipoint_globe)=vsf_in(ix,jy)
               ENDIF
            ENDIF
            ice(ipoint_globe)=ice_in(ix,jy,1)
c
c     If there are no data available for the depth of ice, then
c     assume a depth of 2 metres.  This is based on the default
c     behavior in HadGEM3-A when running with AMIP-2 sea ice.
c
c     NPK 16/12/09 - R3
c            
            IF (.NOT. L_CLIM_ICE_DEPTH) THEN
               IF (ice(ipoint_globe) .GE. 0) THEN
                  icedepth(ipoint_globe)=2.00
               ELSE
                  icedepth(ipoint_globe)=0.00
               ENDIF
            ELSE
               icedepth(ipoint_globe)=icedepth_in(ix,jy,1)
            ENDIF
c
c     If there are no data available for the amount of snow on the ice,
c     then assume that there isn't any.  This is, again, based on the
c     default behavior in HadGEM3-A when running with AMIP-2 sea ice.
c
c     NPK 16/12/09 - R3
c
            IF (.NOT. L_CLIM_SNOW_ON_ICE) THEN
               snowdepth(ipoint_globe)=0.00
            ELSE
               snowdepth(ipoint_globe)=snowdepth_in(ix,jy,1)
            ENDIF
         ENDDO
      ENDDO
      IF (L_OUTKELVIN) SST = SST+TK0
      WRITE(il_mparout,*) 'KPP: Finished creating coupling outputs'
      WRITE(nuout,*) 'KPP: Finished creating coupling outputs'

      DO i=1,jpfldout
c     
c     Export each field to OASIS.  Use a SELECT CASE block to avoid
c     repeated bits of code.  Use the "temporary" variable to transfer
c     the SST and ICE fields to the OASIS "ip_realwp_p" TYPE.
c     
         SELECT CASE (cl_writ(i))
         CASE('OCN_SST')
#ifdef TOYCLIM
            temporary=SST
#else
            CALL ONED_GLOBAL_TWOD_GLOBAL(SST,temporary)
#endif
         CASE('OFRZN01')
#ifdef TOYCLIM
            temporary=ICE
#else
            CALL ONED_GLOBAL_TWOD_GLOBAL(ICE,temporary)
#endif
         CASE ('OSNWTN01')
#ifdef TOYCLIM
            temporary=SNOWDEPTH
#else
            CALL ONED_GLOBAL_TWOD_GLOBAL(SNOWDEPTH,temporary)
#endif            
         CASE ('OHICN01')
#ifdef TOYCLIM
            temporary=ICEDEPTH
#else
            CALL ONED_GLOBAL_TWOD_GLOBAL(ICEDEPTH,temporary)
#endif
         CASE ('SUNOCEAN')
#ifdef TOYCLIM
            temporary=SURF_CURR_X
#else
            IF (L_COUPLE_CURRENTS) THEN
               CALL ONED_GLOBAL_TWOD_GLOBAL(SURF_CURR_X,temporary)
            ELSE
               temporary=0.
            ENDIF
#endif
         CASE ('SVNOCEAN')
#ifdef TOYCLIM
            temporary=SURF_CURR_Y
#else
            IF (L_COUPLE_CURRENTS) THEN
               CALL ONED_GLOBAL_TWOD_GLOBAL(SURF_CURR_Y,temporary)
            ELSE
               temporary=0.
            ENDIF
#endif
         CASE DEFAULT
c
c     The user should never see this - if they do, it means that the model
c     is exporting a field that is not defined here (and likely not defined
c     in the "namcouple" file either).
c
            WRITE(nuout,*) 'Unexpected CASE DEFAULT for i=',i
         END SELECT
         WRITE(nuout,*) 'KPP: Calling PRISM_Put_Proto ',
     +        'for variable ',cl_writ(i)         
         CALL prism_put_proto(il_var_id_out(i),
     +        time_in_seconds,temporary,ierror)
         IF (ierror.NE.PRISM_Ok.and.ierror.LT.PRISM_Sent) THEN
            WRITE(nuout,*) 'KPP: Received error from ',
     +           'PRISM_Put_Proto =',ierror,' sending variable ',
     +           cl_writ(i),' at model time ',time_in_seconds,' sec.'
            WRITE(nuout,*) 'KPP: Aborting coupled integration ...'
            CALL prism_abort_proto(il_comp_id,'couple_io_oasis3.f',
     +           'send')
         ELSE
            WRITE(nuout,*) 'KPP: Successfully called ',
     +           'PRISM_Put_Proto for variable ',cl_writ(i)
         ENDIF
      ENDDO
c
c     End
c
      RETURN
      END SUBROUTINE mpi1_oasis3_output

      SUBROUTINE mpi1_oasis3_terminate
c     
c     Allows KPP to terminate its role in the coupled integration
c     Essentially a wrapper around <prism_terminate_proto>
c
c     NPK 28/09/09 - R3
c     
c     Use OASIS3 (PRISM) modules - the compiler must be able to find
c     these via the Makefile.
c      
      USE mod_kinds_model 
      USE mod_prism_proto
c
      IMPLICIT NONE
c
#include <kpp_oasis3.inc>
c
c     Local variables
c
      INTEGER ierror
c
c     Terminate the integration
c
      WRITE(6,*) 'KPP : Calling prism_terminate_proto(ierror)'
      CALL prism_terminate_proto(ierror)
      WRITE(6,*) 'KPP : Called prism_terminate_proto(ierror)'
      IF (ierror .NE. PRISM_Ok) THEN
         WRITE(il_mparout,*) 'KPP: Received error from ',
     +        'PRISM_Terminate_Proto =',ierror,
     +        ' when terminating model.'
      ENDIF
      WRITE(il_mparout,*) '---- End of the KPP integration ----'
      CLOSE(il_mparout)
c     
c     End
c
      RETURN
      END SUBROUTINE mpi1_oasis3_terminate

      SUBROUTINE mpi1_oasis3_read_netcdf(solar,non_solar,PminusE,
     +     ustress,vstress)
      IMPLICIT NONE

#include <parameter.inc>
#include <netcdf.inc>
     
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
      REAL solar(NPTS),non_solar(NPTS),PminusE(NPTS),
     +     ustress(NPTS),vstress(NPTS)
      REAL*4 temp_in(NPTS)
      INTEGER ncid,status,dimid,varids(6)

      status=NF_OPEN('kpp_atmos_fields.nc',NF_NOWRITE,ncid)
      IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'KPP : Reading coupling fields from input file'
      
      status=NF_INQ_VARID(ncid,'solar',varids(1))
      status=NF_INQ_VARID(ncid,'non_solar',varids(2))
      status=NF_INQ_VARID(ncid,'PminusE',varids(3))
      status=NF_INQ_VARID(ncid,'ustress',varids(4))
      status=NF_INQ_VARID(ncid,'vstress',varids(5))

      status=NF_GET_VARA_REAL(ncid,varids(1),(/1/),(/NPTS/),
     +     temp_in)
      solar=temp_in
      status=NF_GET_VARA_REAL(ncid,varids(2),(/1/),(/NPTS/),
     +     temp_in)
      non_solar=temp_in
      status=NF_GET_VARA_REAL(ncid,varids(3),(/1/),(/NPTS/),
     +     temp_in)
      PminusE=temp_in
      status=NF_GET_VARA_REAL(ncid,varids(4),(/1/),(/NPTS/),
     +     temp_in)
      ustress=temp_in
      status=NF_GET_VARA_REAL(ncid,varids(5),(/1/),(/NPTS/),
     +     temp_in)
      vstress=temp_in

      status=NF_CLOSE(ncid)

      RETURN
      END SUBROUTINE mpi1_oasis3_read_netcdf      

      SUBROUTINE mpi1_oasis3_write_netcdf(solar,non_solar,PminusE,
     +     ustress,vstress)
      IMPLICIT NONE

#include <parameter.inc>
#include <netcdf.inc>
      
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
      REAL solar(NPTS),non_solar(NPTS),PminusE(NPTS),
     +     ustress(NPTS),vstress(NPTS)
      REAL*4 temp_out(NPTS)
      INTEGER ncid,status,dimid,varids(6)

      status=NF_CREATE('kpp_atmos_fields.nc',NF_CLOBBER,ncid)
      IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'KPP : Coupling fields output file created'
      
      status=NF_DEF_DIM(ncid,'point',NPTS,dimid)
      status=NF_DEF_VAR(ncid,'solar',NF_FLOAT,1,(/dimid/),varids(1))
      status=NF_DEF_VAR(ncid,'non_solar',NF_FLOAT,1,(/dimid/),
     +     varids(2))
      status=NF_DEF_VAR(ncid,'PminusE',NF_FLOAT,1,(/dimid/),varids(3))
      status=NF_DEF_VAR(ncid,'ustress',NF_FLOAT,1,(/dimid/),varids(4))
      status=NF_DEF_VAR(ncid,'vstress',NF_FLOAT,1,(/dimid/),varids(5))
      
      status=NF_ENDDEF(ncid)
      
      temp_out=solar
      status=NF_PUT_VARA_REAL(ncid,varids(1),temp_out)
      temp_out=non_solar
      status=NF_PUT_VARA_REAL(ncid,varids(2),temp_out)
      temp_out=PminusE
      status=NF_PUT_VARA_REAL(ncid,varids(3),temp_out)
      temp_out=ustress
      status=NF_PUT_VARA_REAL(ncid,varids(4),temp_out)
      temp_out=vstress
      status=NF_PUT_VARA_REAL(ncid,varids(5),temp_out)
      
      status=NF_CLOSE(ncid)

      RETURN
      END SUBROUTINE mpi1_oasis3_write_netcdf
            
      SUBROUTINE ONED_GLOBAL_ONED_REGIONAL(global,regional)
c
c     Transforms a one-dimensional global field to a one-dimensional
c     regional field.
c     NPK 19/9/09 - R3
c     
      IMPLICIT NONE

#include <parameter.inc>
#include <couple.com>
      
      REAL global(NPTS_GLOBE),regional(NPTS)
      INTEGER*4 ix,jy,ipoint,ipoint_globe

      DO jy=jfirst,jlast
         DO ix=ifirst,ilast
            ipoint_globe=jy*NX_GLOBE+ix
            ipoint=(jy-jfirst)*nx+(ix-ifirst)+1
            regional(ipoint)=global(ipoint_globe)
         ENDDO
      ENDDO
         
      RETURN
      END

      SUBROUTINE TWOD_GLOBAL_ONED_REGIONAL(global,regional)
c
c     Transforms a two-dimensional global field (e.g., one received
c     from HadGEM3-A) to a one-dimensional regional field (as required
c     for KPP).
c     NPK 28/9/09 - R3
c
      IMPLICIT NONE
#include <parameter.inc>
#include <couple.com>
      
      REAL global(NX_GLOBE,NY_GLOBE),regional(NPTS)
      INTEGER*4 ix,jy,ipoint,ipoint_globe

      DO jy=jfirst,jlast
         DO ix=ifirst,ilast
            ipoint=(jy-jfirst)*NX+(ix-ifirst)+1
            regional(ipoint)=global(ix,jy)
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE ONED_GLOBAL_TWOD_GLOBAL(oned,twod)
      IMPLICIT NONE
c
c     Transforms a one-dimensional global field (as produced by
c     the combination of a KPP regional field and a global climatology)
c     to a two-dimensional global field (as required by HadGEM3-A).
c     NPK 28/9/09 - R3
c
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
      END      

#endif /*OASIS3*/
#endif /*COUPLE*/
