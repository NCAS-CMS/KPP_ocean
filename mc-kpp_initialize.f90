  SUBROUTINE MCKPP_INITIALIZE(kpp_3d_fields,kpp_const_fields)
    
    IMPLICIT NONE

    ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
#include <landsea.com>
#include <constants.com>
#include <times.com>
#include <timocn.com>
#include <location.com>
#include <vert_pgrid.com>
#include <proc_swit.com>
#include <proc_pars.com>
#include <initialcon.com>
#include <ocn_advec.com>
#include <ocn_state.com>
#include <ocn_paras.com>
#include <ice_paras.com>
#include <flx_paras.com>
#include <flx_in.com>
#include <output.com>
#include <couple.com>
#include <sstclim.com>
#include <fcorr_in.com>
#include <sfcorr_in.com>
#include <relax_3d.com>
#include <bottomclim.com>
#include <currclim.com>

    !     Input/output
    TYPE(kpp_3d_type),intent(inout) :: kpp_3d_fields
    TYPE(kpp_const_type),intent(inout) :: kpp_const_fields

    ! Local variables    
    INTEGER nuout,nuerr
    PARAMETER(nuout=6,nuerr=0)
    CHARACTER*50 :: output_file,mean_output_file,min_output_file,max_output_file    
    REAL :: alat,alon,delta_lat,delta_lon,dscale
    INTEGER :: i,j,k,l,ipt,ix,iy
    CHARACTER*50 :: forcing_file          

    NAMELIST/NAME_CONSTANTS/grav,vonk,sbc,twopi,onepi,TK0,spd,dpy,&
         epsw,albocn,EL,SL,FL,FLSN
    NAMELIST/NAME_PROCSWIT/LKPP,LRI,LDD,LICE,&
         LBIO,LNBFLX,LTGRID,LRHS,L_SSref
    NAMELIST/NAME_DOMAIN/DMAX,alon,alat,delta_lat,delta_lon,&
         L_STRETCHGRID,dscale,L_REGGRID,L_VGRID_FILE,vgrid_file
    NAMELIST/NAME_START/ L_INITDATA,initdata_file,L_INTERPINIT,&
         L_RESTART,restart_infile
    NAMELIST/NAME_TIMES/ dtsec,startt,finalt,ndtocn
    NAMELIST/NAME_ADVEC/ L_ADVECT,advect_file,L_RELAX_SST,&
         relax_sst_in,relax_sal_in,L_RELAX_CALCONLY,L_RELAX_SAL,&
         L_RELAX_OCNT,relax_ocnt_in
    NAMELIST/NAME_PARAS/ paras_file,L_JERLOV
    NAMELIST/NAME_OUTPUT/ ndt_varout_inst,ndt_singout_inst,output_file,&
         L_RESTARTW,restart_outfile,ndt_varout_mean,ndt_singout_mean,&
         L_OUTPUT_MEAN,L_OUTPUT_INST,L_OUTPUT_RANGE,ndt_per_file,&
         ndt_per_restart,ndt_varout_range,ndt_singout_range,zprof_varout_inst,&
         zprof_varout_mean,zprof_varout_range,zprofs
    NAMELIST/NAME_FORCING/ L_FLUXDATA,forcing_file,L_FCORR_WITHZ,&
         fcorrin_file,ndtupdfcorr,L_VARY_BOTTOM_TEMP,ndtupdbottom,&
         bottomin_file,L_FCORR,L_UPD_FCORR,L_UPD_BOTTOM_TEMP,L_REST,&
         L_PERIODIC_FCORR,L_PERIODIC_BOTTOM_TEMP,fcorr_period,L_SFCORR_WITHZ,&
         sfcorrin_file,ndtupdsfcorr,L_SFCORR,L_UPD_SFCORR,L_PERIODIC_SFCORR,&
         sfcorr_period,bottom_temp_period,sal_file,L_UPD_SAL,L_PERIODIC_SAL,&
         sal_period,ndtupdsal,ocnt_file,L_UPD_OCNT,L_PERIODIC_OCNT,ocnt_period,&
         ndtupdocnt,L_NO_FREEZE,L_NO_ISOTHERM,isotherm_bottom,isotherm_threshold,&
         L_DAMP_CURR,dtuvdamp,L_INTERP_OCNT,ndt_interp_ocnt,L_INTERP_SAL,ndt_interp_sal
    NAMELIST/NAME_COUPLE/ L_COUPLE,ifirst,ilast,jfirst,jlast,L_CLIMSST,sstin_file,&
         L_UPD_CLIMSST,ndtupdsst,L_CPLWGHT,cplwght_file,icein_file,L_CLIMICE,L_UPD_CLIMICE,&
         ndtupdice,L_CLIM_ICE_DEPTH,L_CLIM_SNOW_ON_ICE,L_OUTKELVIN,L_COUPLE_CURRENTS,&
         currin_file,L_CLIMCURR,L_UPD_CLIMCURR,ndtupdcurr,L_PERIODIC_CLIMICE,L_PERIODIC_CLIMSST,&
         climsst_period,climice_period
    NAMELIST/NAME_LANDSEA/ L_LANDSEA,landsea_file
    
    ! This is a bug fix for the IBM XLF compiler, which otherwise complains
    ! about "incorrect characters" in the namelist.  If you are using the
    ! IBM compiler, you need to pass -WF,-DXLF_OLDNAME when you compile MC-KPP.
    ! NPK 5/6/09 - R2
#ifdef XLF_OLDNAME
    CALL SETRTEOPTS("namelist=old")
#endif
        
    allocate(kpp_const_fields%wmt(0:891,0:49))
    allocate(kpp_const_fields%wst(0:891,0:49))
    allocate(kpp_const_fields%tri(0:NZtmax,0:1,NGRID))
    
    ! Open the namelist
    OPEN(75,FILE='3D_ocn.nml')
    
    ! Initialse and read the constants name list
    spd=86400.                ! secs/day
    dpy=360.                  ! days/year
    twopi=8*atan(1.)          ! 2pi
    onepi=twopi/2.            ! pi
    grav=9.816                ! gravity
    vonk=0.4                  ! Von Karman's constant
    TK0=273.15                ! Kelvin of 0degC
    sbc=5.67e-8               ! Stefan Boltzmann Constant
    epsw=1.0                  ! cor.fac for departure of H2O from B.body
    albocn=0.06               ! albedo for seawater
    sice=4.0                  ! salinity of ice(?)
    EL=2.50e6                 ! Latent heat of evap. at 0C (or constant)
    SL=2512200.               ! Latent heat of evap for ice
    FL=334000.                ! Latent heat of fusion for ice
    FLSN=FL                   ! Latent heat of fusion for snow
    READ(75,NAME_CONSTANTS)
    WRITE(nuout,*) 'KPP : Read Namelist CONSTANTS'
!
    ! Initialize and read the processes namelist
    LKPP=.TRUE.
    LRI=.TRUE.
    LDD=.FALSE.
    LICE=.FALSE.
    LBIO=.FALSE.
    LTGRID=.FALSE.
    LNBFLX=.FALSE.
    LRHS=.FALSE.
    L_SSref=.TRUE.
    READ(75,NAME_PROCSWIT)
    WRITE(nuout,*) 'KPP : Read Namelist PROCSWIT'

    ! Initilalize and read the location name list
    DMAX=0.0
    alat=0.0
    alon=0.0
    delta_lat=2.5
    delta_lon=3.75
    dscale=0.0
    L_STRETCHGRID=.FALSE.
    L_REGGRID=.TRUE.
    L_VGRID_FILE=.FALSE.
    READ(75,NAME_DOMAIN)
    IF (DMAX .LE. 0.0) THEN 
       WRITE(nuerr,*) 'KPP : You must specify a depth for the domain'
       CALL MCKPP_ABORT
    ENDIF
    IF ((L_STRETCHGRID) .AND. (dscale .EQ. 0.0)) THEN
       WRITE(nuerr,*) 'KPP : You cannot have dscale=0 for stretched ',&
            'grids'
       CALL MCKPP_ABORT
    ENDIF
    write(nuout,*) 'KPP : Read Namelist DOMAIN'

    ! Initialize and read the landsea name list
    L_LANDSEA=.FALSE.
    READ(75,NAME_LANDSEA)
    WRITE(nuout,*) 'KPP : Read Namelist LANDSEA'
    IF (L_LANDSEA) THEN
       kpp_3d_fields%dlat(1)=alat
       kpp_3d_fields%dlon(1)=alon
       CALL MCKPP_INITIALIZE_LANDSEA(kpp_3d_fields)
    ELSEIF (L_REGGRID) THEN
       DO iy=1,ny
          DO ix=1,nx
             ipt=(iy-1)*nx+ix
             kpp_3d_fields%dlat(ipt)=alat+(iy-1)*delta_lat
             kpp_3d_fields%dlon(ipt)=alon+(ix-1)*delta_lon
             kpp_3d_fields%ocdepth(ipt)=-10000.
             kpp_3d_fields%L_OCEAN(ipt)=.TRUE.
          ENDDO
       ENDDO
    ELSEIF (.NOT. L_REGGRID .AND. .NOT. L_LANDSEA) THEN
       WRITE(nuerr,*) 'KPP : If you set L_REGGRID=.FALSE., you must',&
            ' specify a land-sea mask file from which to read',&
            ' the locations of the gridpoints in the horizontal.'
    ENDIF

    ! If coupling to the GFS, also read in the global land/sea mask on the GFS grid.
    ! This allows KPP to get the global latitudes and longitudes, which it needs
    ! to create a global GRIB file of SSTs to give back to the GFS.
    ! NPK June 2009 - R2
    
#ifdef COUPLE
#ifdef CFS
    IF (L_LANDSEA) CALL read_landsea_global
#endif
#endif
    
    ! Initialize the vertical grid
    CALL MCKPP_INITIALIZE_GEOGRAPHY(L_STRETCHGRID,dscale,kpp_3d_fields,kpp_const_fields)
    WRITE(6,*) 'after init_env zm = ',kpp_const_fields%zm
    
    ! Initialize and read the start name list
    L_INITDATA= .TRUE.
    L_INTERPINIT= .TRUE.
    L_RESTART= .FALSE.
    WRITE(restart_infile,*) 'fort.30'
    READ(75,NAME_START) 
    write(nuout,*) 'KPP : Read Namelist START'
   
    ! Initialize and read the times namelist
    ndtocn=1
    dtsec=0.0
    startt=-999.999
    finalt=-999.999
    READ(75,NAME_TIMES) 
    IF ((dtsec .LE. 0.0) .OR. (startt .LT. 0.0) .OR. (finalt .LT. 0.0)) THEN 
       WRITE(nuerr,*) 'KPP : You must specify values of dtsec,startt,finalt in the namelist'
       CALL MCKPP_ABORT
    ENDIF
    kpp_const_fields%ndtocn=ndtocn
    kpp_const_fields%spd=spd
    kpp_const_fields%dtsec=dtsec
    kpp_const_fields%startt=startt*kpp_const_fields%spd
    kpp_const_fields%finalt=finalt*kpp_const_fields%spd
    kpp_const_fields%dto=kpp_const_fields%dtsec/float(kpp_const_fields%ndtocn)
    kpp_const_fields%nend=int((kpp_const_fields%finalt-kpp_const_fields%startt)/kpp_const_fields%dtsec)
    kpp_const_fields%nstart=nint(kpp_const_fields%startt)/kpp_const_fields%dto
    IF (float(kpp_const_fields%nend*kpp_const_fields%ndtocn) .NE. &
         (kpp_const_fields%finalt-kpp_const_fields%startt)/kpp_const_fields%dto) THEN
       WRITE(nuerr,*) 'KPP : The integration length is not a multiple of the ocean timestep' 
       WRITE(nuerr,*) 'dto=',kpp_const_fields%dto
       WRITE(nuerr,*) 'finalt=',kpp_const_fields%finalt
       WRITE(nuerr,*) 'startt=',kpp_const_fields%startt
       CALL MCKPP_ABORT
    ENDIF
    kpp_const_fields%startt=kpp_const_fields%startt/kpp_const_fields%spd
    kpp_const_fields%finalt=kpp_const_fields%finalt/kpp_const_fields%spd
    kpp_const_fields%time=kpp_const_fields%startt
    WRITE(nuout,*) 'KPP : Read Namelist TIMES'
   
    ! Initialize and read the couple namelist
#ifdef COUPLE
    L_COUPLE=.TRUE.
#else
    L_COUPLE=.FALSE.
#endif
    L_COUPLE_CURRENTS=.FALSE.
    L_OUTKELVIN=.FALSE.
    L_UPD_CLIMSST=.FALSE.
    L_UPD_CLIMICE=.FALSE.
    L_CLIMICE=.FALSE.
    L_CLIMSST=.FALSE.
    L_CLIMCURR=.FALSE. 
    L_BAD_ICE_DEPTH=.FALSE.
    ifirst=1
    ilast=nx
    jfirst=1
    jfirst=ny
    READ(75,NAME_COUPLE)
    write(nuout,*) 'KPP : Read Namelist COUPLE'
    
    ! If the model is coupled or if coupling weights
    ! have been explicitly enabled, initialize the weights.
    ! NPK 10/9/07 - R1
    ! NPK 2/11/09 - Added #ifdef - R3
   
#ifdef COUPLE
    CALL MCKPP_INITIALIZE_COUPLINGWEIGHT(kpp_3d_fields)
#else
    IF (L_CPLWGHT) CALL MCKPP_INITIALIZE_COUPLINGWEIGHT(kpp_3d_fields)
#endif
    ! Initialize and read the advection namelist
    L_ADVECT=.FALSE.
    L_RELAX_SST=.FALSE.
    L_RELAX_CALCONLY=.FALSE.
    DO iy=1,ny
       relax_sst_in(iy)=0.0
       relax_sal_in(iy)=0.0
    ENDDO
    READ(75,NAME_ADVEC)
    IF (L_ADVECT) THEN
       CALL MCKPP_INITIALIZE_ADVECTION(kpp_3d_fields)
    ELSE
       DO ipt=1,npts
          kpp_3d_fields%nmodeadv(ipt,1)=0
          kpp_3d_fields%nmodeadv(ipt,2)=0
       ENDDO
       write(nuout,*) 'KPP : No advection has been specified'
    ENDIF
    write(nuout,*) 'KPP : Read Namelist ADVEC'
    IF (L_RELAX_SST .OR. L_RELAX_SAL .OR. L_RELAX_OCNT) THEN
       CALL MCKPP_INITIALIZE_RELAXATION(kpp_3d_fields,kpp_const_fields)
    ENDIF
    ! Initialize and read the paras namelist
    paras_file='3D_ocnparas.nc'
    L_JERLOV=.TRUE.
    READ(75,NAME_PARAS)
    CALL MCKPP_INITIALIZE_OPTICS(kpp_3d_fields)
    write(nuout,*) 'KPP : Read Namelist PARAS'
    
    ! Initialize and read the forcing namelist
    L_FLUXDATA=.FALSE.
    L_FCORR_WITHZ=.FALSE.
    L_FCORR=.FALSE.
    L_UPD_FCORR=.FALSE.
    L_SFCORR_WITHZ=.FALSE.
    L_SFCORR=.FALSE.
    L_UPD_SFCORR=.FALSE.
    L_UPD_SAL=.FALSE.
    L_VARY_BOTTOM_TEMP=.FALSE.
    L_UPD_BOTTOM_TEMP=.FALSE.
    L_REST=.FALSE.
    L_NO_FREEZE=.FALSE.
    L_NO_ISOTHERM=.FALSE.
    L_DAMP_CURR=.FALSE.
    forcing_file='1D_ocean_forcing.nc'
    ocnT_file='none'
    READ(75,NAME_FORCING)
    write(nuout,*) 'KPP : Read Namelist FORCING'    
    IF (L_FCORR_WITHZ .AND. L_FCORR) THEN
       WRITE(nuerr,*) 'KPP : L_FCORR and L_FCORR_WITHZ are '&
            //'mutually exclusive.  Choose one or neither.'
       CALL MCKPP_ABORT
    ENDIF
    IF (L_SFCORR_WITHZ .AND. L_SFCORR) THEN
       WRITE(nuerr,*) 'KPP : L_SFCORR and L_SFCORR_WITHZ are '&
            //'mutually exclusive.  Choose one or neither.'
       CALL MCKPP_ABORT
    ENDIF
    IF (L_FCORR_WITHZ .AND. L_RELAX_SST) THEN
       WRITE(nuerr,*) 'KPP : L_FCORR_WITHZ and L_RELAX_SST are '&
            //'mutually exclusive.  Choose one or neither.'
       CALL MCKPP_ABORT
    ENDIF
    IF (L_NO_ISOTHERM .AND. (ocnT_file .eq. 'none' .or.&
         sal_file .eq. 'none')) THEN
       WRITE(nuerr,*) 'KPP : If you specify L_NO_ISOTHERM for '&
            //'reseting of isothermal points, you must specify files '&
            //'from which to read climatological ocean temperature '&
            //'(ocnT_file) and salinity (sal_file).'
       CALL MCKPP_ABORT
    ELSEIF (L_NO_ISOTHERM) THEN
       kpp_const_fields%iso_bot=isotherm_bottom
       kpp_const_fields%iso_thresh=isotherm_threshold
    ENDIF
    IF (L_DAMP_CURR) THEN
       kpp_const_fields%dt_uvdamp=dtuvdamp
    ENDIF
            
    ! Initialize and read the output name list
    ndt_varout_inst(:)=0
    ndt_varout_mean(:)=0
    ndt_varout_range(:)=0
    zprof_varout_inst(:)=0
    zprof_varout_mean(:)=0
    zprof_varout_range(:)=0
    ndt_singout_inst(:)=0
    ndt_singout_mean(:)=0
    ndt_singout_range(:)=0     
    zprofs(:,:)=0
    
    L_OUTPUT_MEAN=.FALSE.
    L_OUTPUT_INST=.TRUE.
    L_RESTARTW=.TRUE.      
    kpp_const_fields%ndt_per_restart=kpp_const_fields%nend*kpp_const_fields%ndtocn    
    output_file='KPPocean'
    mean_output_file='KPPocean'
    min_output_file='KPPocean'
    max_output_file='KPPocean'
    
    ! Set up defaults for ndt_per_file (timesteps between creating
    ! new output files) depending on whether and how KPP is coupled.
    ! GFS defaults to one day because the "coupler" operates via
    ! GRIB files that need to be written/read at each coupling timestep.
    ! NPK June 2009 - R2
#ifdef COUPLE
#ifdef CFS
    ndt_per_file=INT(kpp_const_fields%spd/(kpp_const_fields%dtsec/kpp_const_fields%ndtocn))
#else
    ndt_per_file=5*INT(kpp_const_fields%spd/(kpp_const_fields%dtsec/kpp_const_fields%ndtocn))
#endif /*CFS*/
#else
    ndt_per_file=5*INT(kpp_const_fields%spd/(kpp_const_fields%dtsec/kpp_const_fields%ndtocn))
#endif /*COUPLE*/
    READ(75,NAME_OUTPUT)
    write(nuout,*) 'Read Namelist OUTPUT'    
    
    ! Call routine to copy constants and logicals needed for ocean
    ! physics into the kpp_const_fields derived type.  Added for 
    ! compatability with OpenMP DEFAULT(private). NPK 8/2/13
    CALL mckpp_initialize_constants(kpp_const_fields)    

!!! BEYOND THIS POINT, ALL VALUES USED IN OTHER PARTS OF THE CODE
!!! MUST BE A MEMBER OF A DERIVED TYPE OR INCLUDED IN A PARAMETER FILE.
    
    IF (L_CLIMSST) THEN 
       WRITE(6,*) 'MCKPP_INITIALIZE: Calling MCKPP_READ_SST'
       CALL MCKPP_READ_SST(kpp_3d_fields,kpp_const_fields)
       WRITE(6,*) 'MCKPP_INITIALIZE: Returned from MCKPP_READ_SST'
    ENDIF
    IF (L_CLIMICE) CALL MCKPP_READ_ICE(kpp_3d_fields,kpp_const_fields)
    !IF (L_CLIMCURR) CALL read_surface_currents(kpp_3d_fields,kpp_const_fields)
    IF (L_FCORR_WITHZ) THEN 
       WRITE(6,*) 'MCKPP_INITIALIZE: Calling MCKPP_READ_FCORR_3D'
       CALL MCKPP_READ_FCORR_3D(kpp_3d_fields,kpp_const_fields)
       WRITE(6,*) 'MCKPP_INITIALIZE: Returned from MCKPP_READ_FCORR_3D'
    ENDIF
    IF (L_FCORR) CALL MCKPP_READ_FCORR_2D(kpp_3d_fields)
    IF (L_SFCORR_WITHZ) CALL MCKPP_READ_SFCORR_3D(kpp_3d_fields,kpp_const_fields)
    IF (L_SFCORR) CALL MCKPP_READ_SFCORR_2D(kpp_3d_fields)
    IF (L_VARY_BOTTOM_TEMP) CALL MCKPP_READ_TEMPERATURES_BOTTOM(kpp_3d_fields,kpp_const_fields) !,bottom_temp)
    IF (L_RESTART) THEN
       CALL MCKPP_RESTART_IO_READ(kpp_3d_fields,kpp_const_fields,restart_infile)
    ELSE
       WRITE (6,*) 'MCKPP_INITIALIZE: Calling MCKPP_INITIALIZE_OCEAN_PROFILES'
       CALL MCKPP_INITIALIZE_OCEAN_PROFILES(kpp_3d_fields,kpp_const_fields)
       WRITE (6,*) 'MCKPP_INITIALIZE: Returned from MCKPP_INITIALIZE_OCEAN_PROFILES'
       IF (L_UPD_BOTTOM_TEMP) CALL MCKPP_PHYSICS_OVERRIDES_BOTTOMTEMP(kpp_3d_fields,kpp_const_fields)
    ENDIF
    CALL MCKPP_INITIALIZE_FLUXES_VARIABLES(kpp_3d_fields)
    IF (L_RELAX_SAL) CALL MCKPP_READ_SALINITY_3D(kpp_3d_fields,kpp_const_fields)
    IF (L_RELAX_OCNT) CALL MCKPP_READ_TEMPERATURES_3D(kpp_3d_fields,kpp_const_fields)
    IF (L_NO_ISOTHERM .AND. .NOT. L_RELAX_SAL .AND. .NOT. L_RELAX_OCNT) THEN
       CALL MCKPP_READ_TEMPERATURES_3D(kpp_3d_fields,kpp_const_fields)
       CALL MCKPP_READ_SALINITY_3D(kpp_3d_fields,kpp_const_fields)
    ENDIF
    ! Currently, L_INTERP_OCNT implies L_PERIODIC_OCNT to deal with times
    ! before the first time in the input file.
    IF (L_INTERP_OCNT) L_PERIODIC_OCNT=.TRUE.
    IF (L_INTERP_SAL) L_PERIODIC_SAL=.TRUE.

    ! We need to initialize the forcing file (for atmospheric fluxes)
    ! only if KPP is not coupled to an atmospheric model.
    ! NPK June 2009 - R2

#ifndef COUPLE
    kpp_const_fields%forcing_file=forcing_file
    IF (L_FLUXDATA) THEN
       CALL MCKPP_INITIALIZE_FLUXES_FILE(kpp_const_fields)       
    ENDIF
#endif
    
    kpp_const_fields%ntout_vec_inst(:)=1
    kpp_const_fields%ntout_sing_inst(:)=1
    kpp_const_fields%ntout_vec_mean(:)=1
    kpp_const_fields%ntout_sing_mean(:)=1
    kpp_const_fields%ntout_vec_range(:)=1
    kpp_const_fields%ntout_sing_range(:)=1

    kpp_const_fields%zprofs_mask(:,0)=.TRUE.
    kpp_const_fields%zprofs_mask(:,1:N_ZPROFS_MAX)=.FALSE.

    kpp_const_fields%zprofs_nvalid(0)=NZP1
    DO i=1,N_ZPROFS_MAX
       j=1
       DO WHILE (kpp_const_fields%zprofs(j,i) .ne. 0 .and. j .le. NZP1)
          kpp_const_fields%zprofs_mask(kpp_const_fields%zprofs(j,i),i)=.TRUE.
          j=j+1
       END DO
       kpp_const_fields%zprofs_nvalid(i)=j-1
    END DO

    ! Set up the first output files for means and instantaneous fields.
    flen=INDEX(output_file,' ')-1
    kpp_const_fields%day_out=int(kpp_const_fields%startt+(kpp_const_fields%dtsec/ndtocn)*&
         kpp_const_fields%ndt_per_file/kpp_const_fields%spd)
    write(output_file(flen+1:flen+1),'(a)') '_'
    write(output_file(flen+2:flen+6),'(i5.5)') kpp_const_fields%day_out
    write(output_file(flen+7:flen+9),'(3A)') '.nc'
    
    kpp_const_fields%dtout=kpp_const_fields%dto/kpp_const_fields%spd
    IF (kpp_const_fields%L_OUTPUT_INST) THEN
       CALL mckpp_output_initialize(output_file,kpp_3d_fields,kpp_const_fields,'inst')
       CALL mckpp_output_open(output_file,kpp_const_fields%ncid_out)
    ENDIF
    
    IF (kpp_const_fields%L_OUTPUT_MEAN) THEN
       allocate(kpp_3d_fields%VEC_mean(NPTS,NZP1,NVEC_MEAN))
       allocate(kpp_3d_fields%SCLR_mean(NPTS,NSCLR_MEAN))
       flen=INDEX(mean_output_file,' ')-1
       write(mean_output_file(flen+1:flen+1),'(a)') '_'
       write(mean_output_file(flen+2:flen+6),'(i5.5)') kpp_const_fields%day_out
       write(mean_output_file(flen+7:flen+15),'(9A)') '_means.nc'         
       WRITE(nuout,*) 'KPP : Calling init_output for '//mean_output_file
       CALL mckpp_output_initialize(mean_output_file,kpp_3d_fields,&
            kpp_const_fields,'mean')
       CALL mckpp_output_open(mean_output_file,kpp_const_fields%mean_ncid_out)
       kpp_3d_fields%VEC_mean(:,:,:) = 0.
       kpp_3d_fields%SCLR_mean(:,:) = 0.
    ENDIF

    IF (kpp_const_fields%L_OUTPUT_RANGE) THEN       
       allocate(kpp_3d_fields%VEC_range(NPTS,NZP1,NVEC_RANGE,2))
       allocate(kpp_3d_fields%SCLR_range(NPTS,NSCLR_RANGE,2))
       flen=INDEX(min_output_file,' ')-1
       write(min_output_file(flen+1:flen+1),'(a)') '_'
       write(min_output_file(flen+2:flen+6),'(i5.5)') kpp_const_fields%day_out
       write(min_output_file(flen+7:flen+13),'(7A)') '_min.nc'         
       WRITE(nuout,*) 'KPP : Calling init_output for '//min_output_file
       CALL mckpp_output_initialize(min_output_file,kpp_3d_fields,&
            kpp_const_fields,'minx')
       write(max_output_file(flen+1:flen+1),'(a)') '_'
       write(max_output_file(flen+2:flen+6),'(i5.5)') kpp_const_fields%day_out
       write(max_output_file(flen+7:flen+13),'(7A)') '_max.nc'    
       WRITE(nuout,*) 'KPP : Calling init_output for '//max_output_file
       CALL mckpp_output_initialize(max_output_file,kpp_3d_fields,&
            kpp_const_fields,'maxx')
       CALL mckpp_output_open(min_output_file,kpp_const_fields%min_ncid_out)
       CALL mckpp_output_open(max_output_file,kpp_const_fields%max_ncid_out)
       kpp_3d_fields%VEC_range(:,:,:,1)=2E20
       kpp_3d_fields%SCLR_range(:,:,1)=2E20
       kpp_3d_fields%VEC_range(:,:,:,2)=-2E20
       kpp_3d_fields%SCLR_range(:,:,2)=-2E20
    ENDIF

    kpp_const_fields%ntime=0
    CALL mckpp_physics_lookup(kpp_const_fields)
    CALL mckpp_initialize_ocean_model(kpp_3d_fields,kpp_const_fields)
    ! Write out the data from the initial condition
    IF ( .NOT. kpp_const_fields%L_RESTART .AND. kpp_const_fields%L_OUTPUT_INST) THEN
       DO l=1,N_VAROUTS
          WRITE(6,*) l,kpp_const_fields%ndt_varout_inst(l)
          IF (kpp_const_fields%ndt_varout_inst(l) .gt. 0) &
               CALL mckpp_output_inst(kpp_3d_fields,kpp_const_fields,l)               
       ENDDO
       DO l=1,N_SINGOUTS
          IF (kpp_const_fields%ndt_singout_inst(l) .gt. 0)&
               CALL mckpp_output_inst(kpp_3d_fields,kpp_const_fields,l+N_VAROUTS)
       ENDDO
    ENDIF
    
    CLOSE(75)
    RETURN
  END SUBROUTINE mckpp_initialize
