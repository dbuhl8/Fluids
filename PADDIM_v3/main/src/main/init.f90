!+------------------------------------------------------------+
!| This subroutine performs all the necessary initializations |
!| like setting initial condition, setting up the FFT, ...    |
!+------------------------------------------------------------+
SUBROUTINE init(u,Temp,Chem,B,t,dt,dt1,dt2,startstep,restarted)
  USE defprecision_module
  USE MPI
  USE parameter_module, ONLY: Lmax,Mmax,Nmax,Nx,Ny,Nz,dt,dt_initial, &
     &                        FFTW_wisdom_file,uwisdom,              &
     &                        allocate_fft_storage_scheme,           &
     &                        init_fft_storage_scheme,use_FFTW_wisdom_file
  USE mpi_transf_module, ONLY: init_transforms
  USE message_passing_module, ONLY: nprocs1,nprocs2
  USE state_module, ONLY: velocity,buoyancy,field, allocate_uTCB,     & !PH
      &                   init_u_phys,init_Temp_phys,init_Chem_phys, &
      &                   init_B_phys, set_initial_condition
  USE pnetCDF_IO_module, ONLY: pn_read_state_from_dump,pn_read_state_from_dump_hydro,pn_read_state_from_simdat
  IMPLICIT NONE
  TYPE(velocity) :: u
  TYPE(field)    :: B !PH
  TYPE(buoyancy) :: Temp,Chem
  CHARACTER      :: restarted
  REAL(kind=kr)  :: t,dt,dt1,dt2
  INTEGER(kind=ki) :: startstep,ierr

! initialize the parallel FFTs  
  PRINT*,Nx,Ny,Nz,restarted
  IF (use_FFTW_wisdom_file) THEN 
     CALL init_transforms(Lmax,Mmax,Nmax,Nx,Ny,Nz,               &
          &               nprocs1,nprocs2,MPI_COMM_WORLD,        &
          &               FFTW_wisdom_file,uwisdom)
  ELSE
     CALL init_transforms(Lmax,Mmax,Nmax,Nx,Ny,Nz,               &
          &               nprocs1,nprocs2,MPI_COMM_WORLD)
  ENDIF
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
! initialize the FFT storage scheme
  CALL allocate_fft_storage_scheme
  CALL init_fft_storage_scheme

! allocate memory to store the system state
  CALL allocate_uTCB(u,Temp,Chem,B) !PH

! set initial state
  
  IF (restarted.EQ."D" .OR. restarted.EQ."d") THEN 
     ! This is the most common case. reads previous DUMP.
     ! Works as long as DUMP is there, and fields you want to use are in the DUMP
     ! Use this if running hydro run from hydro
     CALL pn_read_state_from_dump(u,Temp,Chem,B,startstep,t,dt) ! PH
     CALL compute_phys_space_vars(u,Temp,Chem,B)
  ELSE
     ! initialize the fields in physical space
     ! This is in the rare cases the DUMP file is corrupted.
     ! Reads data from SIMDAT, as long as file is there, and all fields are in the SIMDAT
     IF (restarted.EQ."S".OR.restarted.EQ."s") THEN
       CALL pn_read_state_from_simdat(u,Temp,Chem,B,startstep,t,dt) ! PH
       PRINT*,"DT=",dt
     ELSE
        IF (restarted.EQ."H" .OR. restarted.EQ."h") THEN
          ! This is for the cases where you want to use the DUMP, but restart 
          ! a magnetic run from a hydro run, while initializing the field from source code
          CALL pn_read_state_from_dump_hydro(u,Temp,Chem,startstep,t,dt) ! PH
          CALL compute_phys_space_vars(u,Temp,Chem,B)
          CALL init_B_phys(B)
        ELSE
        CALL init_u_phys(u)
#ifdef MAGNETIC
        CALL init_B_phys(B)
#endif
#ifdef TEMPERATURE_FIELD
        CALL init_Temp_phys(Temp)
#endif
#ifdef CHEMICAL_FIELD
        CALL init_Chem_phys(Chem)
#endif
        startstep=0
        t = 0._kr
        dt =  dt_initial
        ENDIF
        CALL set_initial_condition(u,Temp,Chem,B) ! PH
     ENDIF
  ENDIF
  dt1 = 0._kr  ! use RK2/CN as starting scheme in multistep AB/BDF3 version
  dt2 = 0._kr  ! there is no meaningful size of previous time steps

END SUBROUTINE init
