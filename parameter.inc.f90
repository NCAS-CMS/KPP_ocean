      INTEGER NZ,NZM1,NZP1,NDIM,NX,NY
      INTEGER*8 npts
      PARAMETER(NZ      = 59  , NZM1 = NZ-1 , NZP1 = NZ+1 )
#ifdef CFS
      PARAMETER(NDIM  =  1 , NX = 87,  NY = 34, NPTS = NX * NY )
#else
#ifdef OASIS2
      PARAMETER(NDIM  =  1 , NX = 129, NY = 73, NPTS = NX * NY )
#else
#ifdef OASIS3
      PARAMETER(NDIM = 1, NX = 192, NY = 112, NPTS = NX*NY)
!     PARAMETER(NDIM = 1, NX = 81, NY = 52, NPTS = NX*NY)
#else
      PARAMETER(NDIM  =  1 , NX = 129, NY = 73, NPTS = NX * NY )
!      PARAMETER(NDIM = 1, NX = 10, NY = 10, NPTS = NX*NY)
#endif /*OASIS3*/
#endif /*OASIS2*/
#endif /*CFS*/
      INTEGER NVEL,NSCLR,NVP1,NSP1,NSB
      PARAMETER(NVEL    =  2 , NSCLR= 2    , NVP1 = NVEL+1,NSP1=NSCLR+1)
      PARAMETER(NSB     =  1 ) !NSCLR-2
      INTEGER itermax
      REAL hmixtolfrac
      PARAMETER(hmixtolfrac = 0.1          , itermax = 200 )
!     temporary grid
      INTEGER NGRID,NZL,NZU,NZDIVmax,NZtmax,NZP1tmax,igridmax
      PARAMETER(NGRID   = NZ )
      PARAMETER(NZL     =  1 , NZU  = 2    , NZDIVmax = 8 )
      PARAMETER(NZtmax  = NZ +(NZL+NZU)*(NZDIVmax-1), NZP1tmax=NZtmax+1)
      PARAMETER(igridmax= 5)
! fluxes and forcing
      INTEGER NSFLXS,NJDT,NSFLXSM1,NSFLXSP2,NDHARM
      PARAMETER(NSFLXS  =  9 , NJDT = 1    , NSFLXSM1 = NSFLXS-1 ,&
           NSFLXSP2 = NSFLXS+2 )
      PARAMETER(NDHARM  =  5 ) 
! ocean advection
      INTEGER maxmodeadv
      PARAMETER(maxmodeadv = 6)
! richardson mixing
      INTEGER MR,MRP1
      PARAMETER(MR      =100 , MRP1 = MR+1)
! rad/conv model
!      INTEGER NPLEV,NPSAVE
!      PARAMETER(NPLEV   = 18 )
!      PARAMETER(NPSAVE  = 10 )
! output buffer
!      INTEGER NDOUT,NBUFF
!      PARAMETER(NDOUT   = 10+NZP1)
!      PARAMETER(NBUFF   = NZP1*(NVp1+NSP1) + NZP1*(NVEL+NSCLR)
!     +                  + NDOUT + 3*NZ + 5*NSFLXS )
!    +                  + 2*NPLEV + NPSAVE ) ! to store rad/conv output
!
!  Parameters for regional coupling
       INTEGER NX_GLOBE,NY_GLOBE
#ifdef CFS
       PARAMETER (NX_GLOBE=192, NY_GLOBE=94)       
#else
#ifdef OASIS2
       PARAMETER (NX_GLOBE=288, NY_GLOBE=217)
#else 
#ifdef OASIS3
       PARAMETER(NX_GLOBE=192, NY_GLOBE=145)
!     PARAMETER (NX_GLOBE=182, NY_GLOBE=152)
#else
       PARAMETER (NX_GLOBE=288, NY_GLOBE=217)
#endif /*OASIS3*/
#endif /*OASIS2*/
#endif /*CFS*/
       INTEGER*8 NPTS_GLOBE
       PARAMETER (NPTS_GLOBE= NX_GLOBE*NY_GLOBE)
!-----------------------------------------------------------------------
! Note in this version:
!    -albedo for ocean is set to 0.06, which is used for QSW from fcomp
!        or fread when rad/conv is not running:
!        albocn=0.06                         (in subroutine init cnsts)
!-----------------------------------------------------------------------
!
! main (permanent) grid
!
!     NZ    : number of layers      
!     NZP1  : number of grid points 
!     NDIM  & NX & NY : dimension of the model(not used in this version)
!     NVEL  : number of velocity components, i.e. 2 for U and V
!     NSCLR : number of scalars, i.e. T, S, and additional scalars
!     NSB   : number of biological scalars, used in "biocommon.inc"
!     hmixtolfrac : convergence tolerance for hmix(new)-hmix(old) 
!             iteration in ocnstep: fraction of layer thickness hm(kmix)
!     itermax : maximum number of hmix iterations (on main or temporary
!             grids.
!
! temporary grid
!
!     NGRID : number of grids = permanent grid + temporary grids,
!             if only p-grid used: NGRID = 1, if t-grid used: NGRID = NZ
!     NZL   & NZU : refinement interval: it is defined on permanent-grid
!             from (kmix+NZU) to (kmix-NZL)
!     NZDIVmax: maximum number of fine(temporary) grid intervals
!             per permanent grid interval. It is needed to dimension the
!             necessary arrays; the number to be used is being read in 
!             as NZDIV.
!     NZtmax: maximum number of layers on temporary grid,
!             the number to be used in this run is read in as NZT.
!     igridmax : maximum number of new grid refinements in ocntgrid
!
! fluxes and forcing
! 
!     NSFLXS: number of fluxes: sflux(NSFLXS,5,0:NJDT)
!     NJDT  : number of older flux values used to extrapolate new fluxes
!     NDHARM: maximum number of harmonics used to specify forcing
!
! ocean advection
!     maxmodeadv: maximum number of different modes for advection
! richardson mixing
!
!     MR    : dimension of "F of Ri" in ri_mix routine
!
! rad/conv model
!
!     NPLEV : number of levels in atmospheric model.
!             (Note: NPLEV=plev necessary in "rad.par")
!     NPSAVE: number of extra atmospheric variables to be saved
!
! dimension of output array
!
!     NDOUT : dimension of output array "dout", so that a number of 
!             scalar parameters, as well as one extra profile 
!             (on layer or grid) can be stored for diagnostic purposes.
!     NBUFF : dimension of data array "buffer"
