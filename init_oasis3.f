#ifdef COUPLE
#ifdef OASIS3
      SUBROUTINE mpi1_oasis3_init(kpp_const_fields)
c
c     Initialize the MPI environment for the model.
c     NPK 13/09/09 - Updated 2/10/09 to support UM 7.1
c     which exchanges 2D fields, whereas the toyclim model
c     exchanges 1D fields.
c
#ifdef OASIS3_MCT
      USE mod_prism
#else
      USE mod_kinds_model
      USE mod_prism_proto
      USE mod_prism_def_partition_proto
      USE mod_prism_put_proto
      USE mod_prism_get_proto
      USE mod_prism_grids_writing
#endif

      IMPLICIT NONE
#include <netcdf.inc>
#include "kpp_3d_type.com"
#include "kpp_oasis3.inc"

      TYPE(kpp_const_type) :: kpp_const_fields
      
      INTEGER nuout,nuerr                 ! KPP output files
      PARAMETER(nuout=6,nuerr=0)
      
      INTEGER il_nbtotproc                ! Total number of processors in MPI communicator
      INTEGER il_commlocal                ! Local communicator for KPP
      INTEGER il_rank                     ! Number of this processor in MPI communicator
      INTEGER il_flag                     ! Flag for writing grids
      INTEGER il_paral(3)                 ! Definitions for mono-process partitions
      INTEGER il_part_id                  ! ID for the local partition (only mono-process supported)
      INTEGER il_var_nodims(2)            ! Rank and "bundles" for the arrays used 
                                          ! for coupling fields (see code below)
#ifdef TOYCLIM /* For the OASIS3 toy model - Exchange 1D fields */
      INTEGER il_var_shape(2)             ! Minimum and maximum array indices for the 
                                          ! coupling fields (see code below)
#else /* For the Unified Model - Exchange 2D fields */
      INTEGER il_var_shape(4)             ! Minimum and maximum array indices for the 
                                          ! coupling fields (see code below)
#endif
      INTEGER i,my_jpfldin

      CHARACTER*6 cp_modnam               ! Name for this model
      CHARACTER*8 choceout                ! Name for the log file for OASIS messages

      INTEGER ierror                      ! Integer error flag
      INTEGER field                       ! Coupling field counter

      ! Initialize coupled environment

      WRITE(cp_modnam,'(A6)') 'toyoce'
      CALL prism_init_comp_proto(il_comp_id, cp_modnam, ierror)
      IF (ierror .NE. PRISM_Ok) THEN
         WRITE(nuout,*) 'KPP: Received error from ',
     +        'PRISM_Init_Comp_Proto = ',ierror
         CALL prism_abort_proto(il_comp_id,'KPP init_oasis3.f','abort1')
         ! Can/should we call MIXED_ABORT here as well?
      ELSE
         WRITE(nuout,*) 'KPP: Successful call to PRISM_Init_Comp_Proto'
      ENDIF

      ! Get local communicator
      CALL prism_get_localcomm_proto(il_commlocal,ierror)
      IF (ierror .NE. PRISM_Ok) THEN
         WRITE(nuout,*) 'KPP: Received error from ',
     +        'PRISM_Get_LocalComm_Proto = ',ierror
      ELSE
         WRITE(nuout,*) 'KPP: Successfully received local communicator'
      ENDIF

      ! Get total number of processors and the number of this processor
      CALL MPI_Comm_Size(il_commlocal,il_nbtotproc,ierror)
      CALL MPI_Comm_Rank(il_commlocal,il_rank,ierror)

      ! Open a log file for OASIS messages
      il_mparout=85+il_rank
      WRITE(choceout,'(A6,I2)') 'KPPout',il_mparout
      OPEN(il_mparout,file=choceout,form='formatted')

      WRITE(il_mparout,*) 'KPP: Number of processes is ',il_nbtotproc
      WRITE(il_mparout,*) 'KPP: Local process number is ',il_rank
      WRITE(il_mparout,*) 'KPP: Local communicator is ',il_commlocal

      ! Define the grids used by KPP (for master processor only)
      IF (il_rank .EQ. 0) THEN
         CALL prism_start_grids_writing(il_flag)
         IF (il_flag .EQ. 1) THEN
      ! Will we ever need to do this?  Do we need to support it?
            WRITE(nuout,*) 'KPP: il_flag=1, so we will write ',
     +           'grids for PRISM'
            CALL prism_terminate_grids_writing()
         ELSE
            WRITE(nuout,*) 'KPP: il_flag/=1, so we will not write ',
     +           'grids for PRISM'
         ENDIF
      ENDIF

      ! Defintions for a mono-process partition
      il_paral ( clim_strategy ) = clim_serial
      il_paral ( clim_offset   ) = 0
      il_paral ( clim_length   ) = NX_GLOBE*NY_GLOBE

      CALL prism_def_partition_proto(il_part_id,il_paral,ierror)
      IF (ierror.NE.PRISM_Ok) THEN
         WRITE(nuout,*) 'KPP: Received error from ',
     +        'PRISM_Def_Partition_Proto = ',ierror
      ELSE
         WRITE(nuout,*) 'KPP: Called PRISM_Def_Partition_Proto'
      ENDIF
      
#ifdef TOYCLIM /* For the OASIS3 toy model - Exchange 1D fields */
      ! Declarations for the coupling field
      il_var_nodims(1) = 1                  ! Rank of the array used for the coupling field
      il_var_nodims(2) = 1                  ! Numbers of the "bundles" (?) in the coupling field (always 1 according to OASIS)
      il_var_shape(1)  = 1                  ! Minimum array index for the coupling field
      il_var_shape(2)  = NX_GLOBE*NY_GLOBE  ! Maximum array index for the coupling field
#else /* For the Unfied Model - Exchange 2D fields */
      il_var_nodims(1) = 2
      il_var_nodims(2) = 1
      il_var_shape(1)  = 1
      il_var_shape(2)  = NX_GLOBE
      il_var_shape(3)  = 1
      il_var_shape(4)  = NY_GLOBE
#endif /*TOYCLIM*/
      
      ! Define the name of each field sent by the KPP model
      ! This needs to be the same name as in the <namcouple> file
      IF (kpp_const_fields%L_COUPLE_FLAGS) THEN
         allocate(kpp_const_fields%cl_writ(
     +        SUM(kpp_const_fields%couple_out_flags)))
         allocate(kpp_const_fields%il_var_id_out(
     +        SUM(kpp_const_fields%couple_out_flags)))
         allocate(kpp_const_fields%cl_read(
     +        SUM(kpp_const_fields%couple_in_flags)))
         allocate(kpp_const_fields%il_var_id_in(
     +        SUM(kpp_const_fields%couple_in_flags)))
      ELSE
         allocate(kpp_const_fields%cl_writ(jpfldout))
         allocate(kpp_const_fields%il_var_id_out(jpfldout))
         allocate(kpp_const_fields%cl_read(jpfldin))
         allocate(kpp_const_fields%il_var_id_in(jpfldin))
         kpp_const_fields%couple_out_flags(:)=1
         kpp_const_fields%couple_in_flags(:)=1
      ENDIF
      kpp_const_fields%il_var_id_in(:)=0
      kpp_const_fields%il_var_id_out(:)=0
      DO i=1,jpfldout           ! Maximum number of possible output fields
         field=1
         IF (kpp_const_fields%couple_out_flags(i) .eq. 1) THEN
            SELECT CASE (i)
            CASE (1)
               kpp_const_fields%cl_writ(kpp_const_fields%fout)='OCN_SST'
            CASE (2)
               kpp_const_fields%cl_writ(kpp_const_fields%fout)='OFRZN01'
            CASE (3)
               kpp_const_fields%cl_writ(kpp_const_fields%fout)=
     +              'OSNWTN01'
            CASE (4)
               kpp_const_fields%cl_writ(kpp_const_fields%fout)='OHICN01'
            CASE (5)
               kpp_const_fields%cl_writ(kpp_const_fields%fout)=
     +              'SUNOCEAN'
            CASE (6)
               kpp_const_fields%cl_writ(kpp_const_fields%fout)=
     +              'SVNOCEAN'
            END SELECT
            kpp_const_fields%fout=kpp_const_fields%fout+1
         ENDIF
      ENDDO
      DO i=1,jpfldin            ! Maximum number of possible input fields
         kpp_const_fields%fin=1
         IF (kpp_const_fields%couple_in_flags(i) .eq. 1) THEN
            SELECT CASE (i)
            CASE(1)
               kpp_const_fields%cl_read(kpp_const_fields%fin)=
     +              'HEATFLUX'
            CASE(2)
#ifdef UM78
               kpp_const_fields%cl_read(kpp_const_fields%fin)='SOLAR'
#endif
#ifdef UM85
               kpp_const_fields%cl_read(kpp_const_fields%fin)='PEN_SOL'
#endif               
            CASE(3)
#ifdef UM78
               kpp_const_fields%cl_read(kpp_const_fields%fin)='WME'
#endif
#ifdef UM85
               IF (kpp_const_fields%L_DIST_RUNOFF) THEN
                  kpp_const_fields%cl_read(kpp_const_fields%fin)=
     +                 'RUNOFF'
               ELSE
                  kpp_const_fields%couple_in_flags(i)=0
               ENDIF
#endif
            CASE(4)
               kpp_const_fields%cl_read(kpp_const_fields%fin)='TRAIN'
            CASE(5)
               kpp_const_fields%cl_read(kpp_const_fields%fin)='TSNOW'
            CASE(6)
               kpp_const_fields%cl_read(kpp_const_fields%fin)='EVAP2D'
            CASE(7)
#ifdef UM78
               kpp_const_fields%cl_read(kpp_const_fields%fin)='LHFLX'
#endif
#ifdef UM85
               kpp_const_fields%cl_read(kpp_const_fields%fin)='TMLT01'
#endif
            CASE(8)
#ifdef UM78
               kpp_const_fields%cl_read(kpp_const_fields%fin)='TMLT01'
#endif
#ifdef UM85
               kpp_const_fields%cl_read(kpp_const_fields%fin)='TAUX'
#endif
            CASE(9)
#ifdef UM78
               kpp_const_fields%cl_read(kpp_const_fields%fin)='BMLT01'
#endif
#ifdef UM85
               kpp_const_fields%cl_read(kpp_const_fields%fin)='TAUY'
#endif
#ifdef UM78
            CASE(10)
               kpp_const_fields%cl_read(kpp_const_fields%fin)='TAUX'
            CASE(11)
               kpp_const_fields%cl_read(kpp_const_fields%fin)='TAUY'
#endif
            END SELECT
            kpp_const_fields%fin=kpp_const_fields%fin+1
         ENDIF
      ENDDO
     

      DO i=1,kpp_const_fields%fout
         CALL prism_def_var_proto(kpp_const_fields%il_var_id_out(i),
     +        kpp_const_fields%cl_writ(i),
     +        il_part_id,il_var_nodims,PRISM_Out,il_var_shape,
     +        PRISM_Real,ierror)
         IF (ierror.NE.PRISM_Ok) THEN
            WRITE(nuout,*) 'KPP: Received error from ',
     +           'PRISM_Def_Var_Proto = ',ierror,'for variable ',
     +           kpp_const_fields%cl_writ(i),' (output field)'
         ELSE
            WRITE(nuout,*) 'KPP: Called PRISM_Def_Var_Proto for ',
     +           'variable ',kpp_const_fields%cl_writ(i),
     +           ' (output field)'
         ENDIF
      ENDDO

      DO i=1,kpp_const_fields%fin
         CALL prism_def_var_proto(kpp_const_fields%il_var_id_in(i),
     +        kpp_const_fields%cl_read(i),
     +        il_part_id,il_var_nodims,PRISM_In,il_var_shape,
     +        PRISM_Real,ierror)
         IF (ierror.NE.PRISM_Ok) THEN
            WRITE(nuout,*) 'KPP: Received error from ',
     +           'PRISM_Def_Var_Proto = ',ierror,'for variable',
     +           kpp_const_fields%cl_read(i),' (input field)'
         ELSE
            WRITE(nuout,*) 'KPP: Called PRISM_Def_Var_Proto for ',
     +           'variable ',kpp_const_fields%cl_read(i),
     +           ' (input field)'
         ENDIF
      ENDDO

      CALL prism_enddef_proto(ierror)
      IF (ierror.NE.PRISM_Ok) THEN
         WRITE(nuout,*) 'KPP: Received error from ',
     +        'PRISM_enddef_proto = ',ierror
      ELSE
         WRITE(nuout,*) 'KPP: Called PRISM_Enddef_Proto'
      ENDIF

      RETURN
            
      END SUBROUTINE mpi1_oasis3_init
#endif /*OASIS3*/
#endif /*COUPLE*/
