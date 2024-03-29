!+--------------------------------------------------------+
!|                      PADDIM-Code                       |
!|                        =====                           |
!|          (Pa)rallel (D)oube (Di)ffusion (MHD) Code     |
!|                                                        |
!| Spectral Code to solve the double diffusive Boussinesq |
!| MHD equations in a triply periodic cube.               |
!+--------------------------------------------------------+
!| Written by:    Stephan Stellmach         (2008-2011)   |
!|                (stellma@uni-muenster.de)               |
!| MHD extension* by:   Peter Harrington    (2018)        |
!|                      (pharring@ucsc.edu)               |
!|    *additions/modifications marked by ! PH comment     |
!| Forcing module* by:  Dante Buhl          (2023)        |
!|                      (dbuhl@ucsc.edu)                  |
!|    *additions/modificaions marked by ! DB comment      |
!|                                                        |
!| To prevent overlap, please always contact me before    |
!| using this code on a particular problem!               |
!| For personal use only. Do not distribute!              |
!+--------------------------------------------------------+
PROGRAM double_diffusion
  USE defprecision_module
  USE message_passing_module, ONLY: start_mpi,stop_mpi,myid
  USE main_module, ONLY: init,read_parameter,free_allocated_memory, &
      &                  timestep_AB_BDF3,timestep_RK2_CN 
  USE state_module, ONLY: velocity, buoyancy, field ! PH
  USE parameter_module, ONLY: nsteps
  USE IO_module, ONLY: open_files,write_compressed_file,write_diagnostics_file,close_files, &
      &                write_output_files
  USE message_passing_module, ONLY : myid
  USE forcing_module, ONLY: read_stochastic_forcing, close_forcing, init_forcing, forcing ! DB
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  CHARACTER :: restarted
  TYPE(velocity) :: u
  TYPE(field)    :: B ! PH
  TYPE(buoyancy) :: Temp,Chem
  INTEGER(kind=ki) :: istep,startstep,error_code = 0
  REAL(kind=kr) :: t,dt,dt1,dt2
  REAL(kind=krd)     :: time_begin=0.,time_end=0.
  
! Initialize MPI
  CALL start_mpi

! Read parameter values for the run
  CALL read_parameter(restarted, error_code) !DB added error_code

  if(error_code .ne. 0) then
    goto 30
  end if

! perform necessary initializations
  CALL init(u,Temp,Chem,B,t,dt,dt1,dt2,startstep,restarted) ! PH

  CALL init_forcing ! DB
 
#ifdef STOCHASTIC_FORCING
  CALL read_stochastic_forcing(restarted, t) ! DB 
#endif
!  write(*,*) 'Done initial conditions'

! open output files 
  CALL open_files

#ifdef STOCHASTIC_FORCING
! Populates forcing arrays before writing diagnostic file on a restart (ensures continuity in power input column)
  if (restarted .ne. "N") then
    CALL forcing(t) !DB
  end if
#endif
!  write(*,*) 'Done opening files'
  CALL write_diagnostics_file(u,Temp,Chem,B,0,t,dt) ! PH
!  write(*,*) 'Done writing diag files'
  CALL write_compressed_file(u,Temp,Chem,0,t,dt)
!  write(*,*) 'Done writing compressed filse'
   
! start time stepping

#ifdef AB_BDF3
  ! Third order Adams-Bashforth / Backward-Differencing multi-step method
  time_stepping_loop : DO istep = startstep+1,startstep+nsteps
     IF (myid.EQ.0) WRITE(*,'(a,I7,a,F8.3,a)') "step :",istep,"  ",time_end-time_begin," s"
     time_begin=MPI_WTIME()
     IF (istep.LT.startstep+3) THEN
        ! use second order Runge-Kutta / Crank Nicholson as starting scheme
        CALL timestep_RK2_CN(u,Temp,Chem,B,istep,t,dt,dt1,dt2) ! PH
     ELSE 
        CALL timestep_AB_BDF3(u,Temp,Chem,B,istep,t,dt,dt1,dt2) ! PH
     ENDIF
     time_end=MPI_WTIME()
     CALL write_output_files(u,Temp,Chem,B,t,dt,istep) ! PH
  END DO time_stepping_loop
#endif

! close output files
  CALL close_files

! closes forcing files, 
  CALL close_forcing ! DB, if stochastic forcing is defined, creates forcing restart files
 !IMPORTANT: if the dump file is not written close to the last timestep, the forcing (and other outputs) may be disconinuous after
 !restart!
 ! This is easy to avoid as long as the number of timesteps is a multiple of restart_info_... in the parameter file. 

! free allocated memory
  CALL free_allocated_memory(u,Temp,Chem,B) ! PH


30     print *, "Error in read_parameter. Check logs to see why" 
        stop

! Stop MPI
  CALL stop_mpi



END PROGRAM double_diffusion
