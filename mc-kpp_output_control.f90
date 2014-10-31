SUBROUTINE mckpp_output_control(kpp_3d_fields,kpp_const_fields,kpp_timer)
  IMPLICIT NONE

#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  TYPE(kpp_timer_type) :: kpp_timer

  CHARACTER(LEN=50) :: output_file,mean_output_file,min_output_file,max_output_file
  INTEGER :: j,flen

  ! Need to handle rotation of all output files

  IF (kpp_const_fields%L_OUTPUT_MEAN) THEN 
     WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Calling MCKPP_OUTPUT_COMPUTE_MEANS'
     !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
     !CALL KPP_TIMER_TIME(kpp_timer,'Time-meaning of output',1)
     CALL MCKPP_OUTPUT_COMPUTE_MEANS(kpp_3d_fields,kpp_const_fields)
     !CALL KPP_TIMER_TIME(kpp_timer,'Time-meaning of output',0)
     !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
     WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Returned from MCKPP_OUTPUT_COMPUTE_MEANS'
  ENDIF
  
  IF (kpp_const_fields%L_OUTPUT_RANGE) THEN
     WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Calling MCKPP_OUTPUT_COMPUTE_RANGES'
     !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
     !CALL KPP_TIMER_TIME(kpp_timer,'Time-meaning of output',1)
     CALL MCKPP_OUTPUT_COMPUTE_RANGES(kpp_3d_fields,kpp_const_fields)
     !CALL KPP_TIMER_TIME(kpp_timer,'Time-meaning of output',0)
     !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
     WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Returned from MCKPP_OUTPUT_COMPUTE_RANGES'
  ENDIF  

  DO j=1,N_VAROUTS
     IF (kpp_const_fields%ndt_varout_mean(j) .gt. 0 .and. kpp_const_fields%L_OUTPUT_MEAN) THEN
        IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndt_varout_mean(j)).EQ.0) THEN 
           !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
           !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
           WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Calling MCKPP_OUTPUT_MEAN'
           CALL MCKPP_OUTPUT_MEAN(kpp_3d_fields,kpp_const_fields,j)
           WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Returned from MCKPP_OUTPUT_MEAN'
           !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
           !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
        ENDIF
     ENDIF
     IF (kpp_const_fields%ndt_varout_inst(j) .gt. 0 .and. kpp_const_fields%L_OUTPUT_INST) THEN
        IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndt_varout_inst(j)).EQ.0) THEN 
           !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
           !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
           WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Calling MCKPP_OUTPUT_INST'
           CALL MCKPP_OUTPUT_INST(kpp_3d_fields,kpp_const_fields,j)
           WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Returned from MCKPP_OUTPUT_INST'
           !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
           !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
        ENDIF
     ENDIF
     IF (kpp_const_fields%ndt_varout_range(j) .gt. 0 .and. kpp_const_fields%L_OUTPUT_RANGE) THEN
        IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndt_varout_range(j)).EQ.0) THEN
           !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
           !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
           WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Calling MCKPP_OUTPUT_RANGE'           
           CALL MCKPP_OUTPUT_RANGE(kpp_3d_fields,kpp_const_fields,j)
           WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Returned from MCKPP_OUTPUT_RANGE'
           !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
           !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
        ENDIF
     ENDIF
  ENDDO
  DO j=1,N_SINGOUTS
     IF (kpp_const_fields%ndt_singout_mean(j) .gt. 0 .and. kpp_const_fields%L_OUTPUT_MEAN) THEN
        IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndt_singout_mean(j)).EQ.0) THEN
           !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
           !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
           WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Calling MCKPP_OUTPUT_MEAN'
           CALL MCKPP_OUTPUT_MEAN(kpp_3d_fields,kpp_const_fields,j+N_VAROUTS)
           WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Returned from MCKPP_OUTPUT_MEAN'          
           !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
           !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
        ENDIF
     ENDIF
     IF (kpp_const_fields%ndt_singout_inst(j) .gt. 0 .and. kpp_const_fields%L_OUTPUT_INST) THEN
        IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndt_singout_inst(j)).EQ.0) THEN 
           !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
           !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
           WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Calling MCKPP_OUTPUT_INST'
           CALL MCKPP_OUTPUT_INST(kpp_3d_fields,kpp_const_fields,j+N_VAROUTS)
           WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Returned from MCKPP_OUTPUT_INST'           
           !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
           !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
        ENDIF
     ENDIF
     IF (kpp_const_fields%ndt_singout_range(j) .gt. 0 .and. kpp_const_fields%L_OUTPUT_RANGE) THEN
        IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndt_singout_range(j)) .EQ. 0) THEN
           !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
           !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
           WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Calling MCKPP_OUTPUT_RANGE'                     
           CALL MCKPP_OUTPUT_RANGE(kpp_3d_fields,kpp_const_fields,j+N_VAROUTS)
           WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: Returned from MCKPP_OUTPUT_RANGE'
           !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
           !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
        ENDIF
     ENDIF
  ENDDO

  !day_out=INT(kpp_const_fields%startt+(kpp_const_fields%dtsec/&
  !     kpp_const_fields%ndtocn)*kpp_const_fields%ndt_per_file/kpp_const_fields%spd)
  WRITE(6,*) 'MCKPP_OUTPUT_CONTROL: day_out=',kpp_const_fields%day_out,' and time=',kpp_const_fields%time
       
  IF (kpp_const_fields%time .GE. REAL(kpp_const_fields%day_out) .AND. kpp_const_fields%ntime .NE. &
       kpp_const_fields%nend*kpp_const_fields%ndtocn) THEN
     kpp_const_fields%day_out=kpp_const_fields%day_out+NINT(kpp_const_fields%dtsec/REAL(kpp_const_fields%ndtocn)*&
          FLOAT(kpp_const_fields%ndt_per_file)/kpp_const_fields%spd)     
     output_file='KPPocean_'
     mean_output_file='KPPocean_'
     min_output_file='KPPocean_'
     max_output_file='KPPocean_'
     IF (kpp_const_fields%L_OUTPUT_INST) THEN
        !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
        !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
        flen=INDEX(output_file,' ')-1
        WRITE(output_file(flen+1:flen+8),'(i5.5,3A)') kpp_const_fields%day_out,'.nc'
        CALL MCKPP_OUTPUT_CLOSE(kpp_const_fields%ncid_out)
        kpp_const_fields%ntout_vec_inst(:)=1
        kpp_const_fields%ntout_sing_inst(:)=1
        CALL MCKPP_OUTPUT_INITIALIZE(output_file,kpp_3d_fields,kpp_const_fields,'inst')
        CALL MCKPP_OUTPUT_OPEN(output_file,kpp_const_fields%ncid_out)
        !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
        !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
     ENDIF
     IF (kpp_const_fields%L_OUTPUT_MEAN) THEN
        !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
        !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)        
        flen=INDEX(mean_output_file,' ')-1
        write(mean_output_file(flen+1:flen+14),'(i5.5,9A)'),kpp_const_fields%day_out,'_means.nc'                 
        CALL MCKPP_OUTPUT_CLOSE(kpp_const_fields%mean_ncid_out)
        kpp_const_fields%ntout_vec_mean(:)=1
        kpp_const_fields%ntout_sing_mean(:)=1
        CALL MCKPP_OUTPUT_INITIALIZE(mean_output_file,kpp_3d_fields,kpp_const_fields,'mean')
        CALL MCKPP_OUTPUT_OPEN(mean_output_file,kpp_const_fields%mean_ncid_out)
        !kpp_const_fields%nout_mean=1
        !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
        !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
     ENDIF
     IF (kpp_const_fields%L_OUTPUT_RANGE) THEN        
        !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
        !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',1)
        flen=INDEX(min_output_file,' ')-1
        write(min_output_file(flen+1:flen+12),'(i5.5,A7)') kpp_const_fields%day_out,'_min.nc'
        flen=INDEX(max_output_file,' ')-1
        write(max_output_file(flen+1:flen+12),'(i5.5,A7)') kpp_const_fields%day_out,'_max.nc'
        CALL MCKPP_OUTPUT_CLOSE(kpp_const_fields%min_ncid_out)
        CALL MCKPP_OUTPUT_CLOSE(kpp_const_fields%max_ncid_out)
        kpp_const_fields%ntout_vec_range(:)=1
        kpp_const_fields%ntout_sing_range(:)=1
        CALL MCKPP_OUTPUT_INITIALIZE(min_output_file,kpp_3d_fields,kpp_const_fields,'minx')
        CALL MCKPP_OUTPUT_INITIALIZE(max_output_file,kpp_3d_fields,kpp_const_fields,'maxx')
        CALL MCKPP_OUTPUT_OPEN(min_output_file,kpp_const_fields%min_ncid_out)
        CALL MCKPP_OUTPUT_OPEN(max_output_file,kpp_const_fields%max_ncid_out)
        !CALL KPP_TIMER_TIME(kpp_timer,'Writing output',0)
        !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
     ENDIF
     STOP
  ENDIF
  

END SUBROUTINE mckpp_output_control

SUBROUTINE mckpp_output_open(file,ncid)      
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
  !#include <kpp_3d_type.com>
  !#include <output.com>

  INTEGER status
  INTEGER,intent(out) :: ncid
  CHARACTER*50,intent(in) :: file      
  
  status=NF_OPEN(file,NF_WRITE,ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  RETURN
END SUBROUTINE mckpp_output_open

SUBROUTINE mckpp_output_close(ncid)      
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
  
#include <netcdf.inc>
  !#include <kpp_3d_type.com>
  !#include <output.com>

  INTEGER status
  INTEGER, intent(in) :: ncid
  
  status=NF_CLOSE(ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  RETURN
END SUBROUTINE mckpp_output_close
