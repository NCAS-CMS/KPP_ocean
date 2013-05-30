      INTEGER N_VAROUTS
      INTEGER N_SINGOUTS
      PARAMETER (N_VAROUTS=23,N_SINGOUTS=10)
      LOGICAL L_VAROUT(N_VAROUTS),L_SINGOUT(N_SINGOUTS)
      LOGICAL L_MEAN_VAROUT(N_VAROUTS),L_MEAN_SINGOUT(N_SINGOUTS)
      LOGICAL L_OUTPUT_MEAN, L_OUTPUT_INST
      LOGICAL L_RESTARTW
      CHARACTER*4 restart_time
      INTEGER ndtout,ndtout_mean,nout,nout_mean,ndt_per_restart
      REAL*4 dtout
      REAL*4 missval
      PARAMETER (missval=1.e20)
************************************************************************
* Common blocks to control outputs, each variable has a number associated 
* with it given below
*
*     VAROUTS
*     1    =    U field     (m/s)
*     2    =    V field     (m/s)
*     3    =    T field     (degC)
*     4    =    S field     (o/oo)
*     5    =    Buoyancy    (m/s2)
*
*     6    =    w'u'        (m2/s2)
*     7    =    w'v'        (m2/s2)
*     8    =    w'T'        (degC m/s)
*     9    =    w'S'        (o/oo m/s)
*     10   =    w'B'        (m2/s3)
*
*     11   =    w'T'(NT)    (degC m/s)
*
*     12   =    difm        (m2/s)
*     13   =    dift        (m2/s)
*     14   =    difs        (m2/s)
*
*     15   =    rho         (kg/m3)
*     16   =    cp          (J/kg/K)
*    
*     17   =    scorr       (o/oo/s)
*     18   =    Rig         (local Richardson number in kpp.f - NPK diagnostic)
*     19   =    dbloc       (local delta buoyancy used to calculate Rig - NPK diagnostic)
*     20   =    shsq        (local shear-squared term used to calculate Rig - NPK diagnostic)
*     21   =    tinc_fcorr  (temperature increment from flux corrections with depth, K)
*     22   =    ocnTcorr    (K/s)
*     23   =    sinc_fcorr  (salinity increment from flux corrections with depth,o/oo)
*
*     SINGOUTS
*     1 = hmix     (m)  : single level field.
*     2 = fcorr    (W/m^2) : single level field
*     3 = taux_in  (N/m^2)
*     4 = tauy_in  (N/m^2)
*     5 = solar_in    (W/m^2)*     6 = nsolar_in   (W/m^2)
*     7 = PminusE_in  (W/m^2)
*     8 = cplwght
*     9 = freeze_flag (unitless) : fraction of levels at which temperature was < -1.8C and was reset to -1.8C.
*    10 = isotherm_flag (unitless)  : 1/0 for whether T/S profile was reset to climatology because of isothermal detection routine 
************************************************************************
      INTEGER ncid_out,mean_ncid_out,day_out,flen,ndt_per_file
      CHARACTER*50 output_file,mean_output_file,restart_outfile
      
      INTEGER londim,latdim,zdim,ddim,hdim,timdim
      INTEGER lon_id,lat_id,z_id,h_id,d_id,time_id

      INTEGER varid(N_VAROUTS),singid(N_SINGOUTS),mean_varid(N_VAROUTS),
     &     mean_singid(N_SINGOUTS)
      INTEGER NVEC_MEAN,NSCLR_MEAN
      PARAMETER(NVEC_MEAN=N_VAROUTS,NSCLR_MEAN=N_SINGOUTS)
            
      common / output / l_varout,l_singout,l_mean_varout,
     &     l_mean_singout,dtout,nout,nout_mean,restart_outfile,
     &     ndtout,ndtout_mean,l_restartw,l_output_mean,l_output_inst,
     &	   ndt_per_restart,restart_time
      common /ncdf_out/ mean_ncid_out,ncid_out,londim,latdim,zdim,ddim,
     &     hdim,timdim,time_id,varid,singid,day_out,output_file,flen,
     &     mean_output_file,mean_varid,mean_singid,ndt_per_file

