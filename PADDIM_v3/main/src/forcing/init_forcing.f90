! File: init_forcing.f90
! Author: Dante Buhl
! Purpose: initializes forcing arrays and reads parameters in from the parameter file.

subroutine init_forcing

    USE defprecision_module
    USE parameter_module
    USE message_passing_module, ONLY: myid
    USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec, &
                                & mysy_phys, myey_phys, mysz_phys, myez_phys
    
#ifdef MPI_MODULE
    USE MPI
#else
    INCLUDE "mpif.h"
#endif

    implicit none
    
    character :: restarted
    integer(kind=ki) :: i, j, ierr
    real(kind=kr) :: hkx, hky
    real(kind=kr) :: max_forced_wavenumber
    !some default parameters
    integer(kind=ki) :: number_of_points_in_window = 10_ki
    integer(kind=ki) :: number_of_points_between_window_update = 250_ki
    real(kind=kr) :: delta_t = 0.5_kr
    real(kind=kr) :: gaussian_timescale = 1._kr
    real(kind=kr) :: relative_eigenvalue_tolerance = 10_kr**-4
    logical :: complexToReal = .true.
    logical :: keepPrevForcing = .false.
   
    NAMELIST /forcing_values/ max_forced_wavenumber
    NAMELIST /forcing_values/ number_of_points_in_window
    NAMELIST /forcing_values/ number_of_points_between_window_update 
    NAMELIST /forcing_values/ delta_t                                
    NAMELIST /forcing_values/ gaussian_timescale                     
    NAMELIST /forcing_values/ relative_eigenvalue_tolerance         
    NAMELIST /forcing_values/ complexToReal
    NAMELIST /forcing_values/ keepPrevForcing

    IF (myid.EQ.0) THEN
    ! Read data from given namelist
    ! -----------------------------
        READ(*,NML=forcing_values,IOSTAT=ierr)
        IF (ierr .NE. 0 ) THEN
            PRINT*,myid,'READ ERROR IN INPUT NAMELIST'
            STOP
        ENDIF
    
        KMAX_forcing                =  max_forced_wavenumber
        window_pts                  =  number_of_points_in_window
        window_skip                 =  number_of_points_between_window_update
        tstep                       =  delta_t
        gaussian_tmscl              =  gaussian_timescale 
        tol                         =  relative_eigenvalue_tolerance
        c2r                         =  complexToReal
        usepf                       =  keepPrevForcing
    
    ENDIF
    ! Broadcast this to the rest of the processors. 
    CALL MPI_BCAST(KMAX_forcing          ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(window_pts            ,1,   MPI_INTEGER   ,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(window_skip           ,1,   MPI_INTEGER   ,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(gaussian_tmscl        ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(tol                   ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(tstep                 ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(c2r                   ,1,   MPI_LOGICAL   ,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(usepf                 ,1,   MPI_LOGICAL   ,0,MPI_COMM_WORLD,ierr)

    allocate(waveNumMap(mysx_spec:myex_spec, mysy_spec:myey_spec))
    ALLOCATE(force_spec(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,vec_x:vec_z))
    ALLOCATE(force_real(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys,vec_x:vec_z))

end subroutine init_forcing
