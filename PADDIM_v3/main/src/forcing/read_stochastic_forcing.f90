! File: read_stochastic_forcing.f90
! Author: Dante Buhl
! Purpose: Read gaussian processes from text file (or generates them ) and place in a
! fortran array

subroutine read_stochastic_forcing(restarted, t)
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
  
  character, intent(in)       :: restarted
  real(kind=kr), intent(in)   :: t
  logical                     :: paramsMatch, forcingexists
  integer(kind=ki)            :: i, j, readNumColumns, ierr, rs, numTrashLines
  real(kind=kr)               :: hkx, hky, khoriz
  real(kind=kr)               :: max_forced_wavenumber, timeStop
  integer(kind=ki)            :: number_of_points_in_window = 10_ki
  integer(kind=ki)            :: number_of_points_between_window_update = 250_ki
  real(kind=kr)               :: delta_t = 0.5_kr
  real(kind=kr)               :: gaussian_timescale = 1._kr
  real(kind=kr)               :: relative_eigenvalue_tolerance = 10_kr**-4
  real(kind=kr)               :: last_timestep_of_forcing
  character*50                :: string

#ifdef STOCHASTIC_FORCING
  c2r = .true.
#endif
 
  headerFormatReal = "('# ', A, F15.8)"
  headerFormatInt = "('# ', A, I6)"


  numGPcolumns = 0_ki
  waveNumMap = 1_ki
  if (n_wrt_dump > Nsteps) then
    numGProws = NINT((Nsteps*dt_max)/tstep, ki) + 20_ki
  else 
    ! sets the length of the forcing arrays equal to lenth between dumps
    numGProws = NINT((n_wrt_dump*dt_max/tstep,ki) + 1_ki
  end if
 
  ! Storing values in waveNumMap
  ! ----------------------------------------------------------------------------------------- 
  if(Mmax .ne. 0) then
    do j=mysy_spec,myey_spec
      hky = ky(j)
      do i=mysx_spec,myex_spec
        hkx = kx(i)
        khoriz = sqrt(hky**2 + hkx**2)
        if((khoriz .le. KMAX_forcing) .and. (hkx .ne. 0.0 .and. hky .ne. 0.0)) then
          numGPcolumns = numGPcolumns + 1_ki
          waveNumMap(i, j) = numGPcolumns
        end if 
      end do
    end do
  else 
    do i = mysx_spec, myex_spec
      if((kx(i) .le. KMAX_forcing) .and. (kx(i) .gt. 0_ki)) then
        numGPcolumns = numGPcolumns + 1_ki
        waveNumMap(i, :) = numGPcolumns
      end if 
    end do
  end if
  ! ----------------------------------------------------------------------------------------- 
  numGPcolumns = numGPcolumns*2
  ! ----------------------------------------------------------------------------------------- 
  forcingCPU = .false.
  ! ----------------------------------------------------------------------------------------- 
  ! need to add an MPI Reduce to get total number of forced wavenumbers
  ! call MPI_REDUCE(tot_num_forced,...,0,MPI_SUM,MPI_INT,ie) 

  ! -----------------------------------------------------------------------------------------
  ! If this processor needs forcing, then do the forcing 
  ! -----------------------------------------------------------------------------------------
  if(numGPcolumns .gt. 0) then
    allocate(gpForcingVals(numGProws, numGPcolumns))
    allocate(gpTimeVals(numGProws)) 
    !allocate(interpSlopeVals(numGProws-1_ki, numGPcolumns))
    dataFormat = "("//trim(str(numGPcolumns+1))//"(F11.5, '    '))"
    forcingCPU = .true.
    timeStop = Nsteps*dt_max + 10.0_kr*dt_max + 20.0_kr*tstep

    ! Make sure there is a forcing file for this processor with enough columns
    ! ------------------------------------------------------------------------------------- 
    !inquire(file="forcing_data/forcing"//trim(str(myid))//".dat", exist=forcingexists, iostat=rs)
    
    ! this entire if statement needs to be tweaked

    !if (.not. forcingexists) then 
      !CALL gaussian(numGPcolumns, window_pts, window_skip, timeStop, tstep, &
        !gaussian_tmscl, tol, myid, numTrashLines, .false.)
      
    !else if (usepf) then
      !!print *, "Using previous forcing"

      !read parameters in from forcing file
      !open(20, file="forcing_data/forcing"//trim(str(myid))//".dat")
        !read(20, headerFormatInt)  string, readNumColumns
        !read(20, headerFormatInt)  string, number_of_points_in_window
        !read(20, headerFormatInt)  string, number_of_points_between_window_update
        !read(20, headerFormatReal) string, gaussian_timescale
        !read(20, headerFormatReal) string, delta_t
        !read(20, headerFormatReal) string, relative_eigenvalue_tolerance
        !read(20, headerFormatReal) string, last_timestep_of_forcing
        !read(20, headerFormatInt)  string, numTrashLines
      !close(20)
      !this set of logic statements might be irrelevant
      !paramsMatch = .true.
      !if (readNumColumns .lt. numGPcolumns) then
        !paramsMatch = .false. 
      !else if (number_of_points_in_window .ne. window_pts) then
        !paramsMatch = .false.
      !else if (number_of_points_between_window_update .ne. window_skip) then
        !paramsMatch = .false.
      !else if (gaussian_timescale .ne. gaussian_tmscl) then
        !paramsMatch = .false.
      !else if (delta_t .ne. tstep) then
        !paramsMatch = .false.
      !else if (relative_eigenvalue_tolerance .ne. tol) then
        !paramsMatch = .false.
      !else if (last_timestep_of_forcing .ne. timeStop) then
        !paramsMatch = .false.
      !end if
      !if (.not. paramsMatch) then
        !print *, "WARNING: Previous Forcing has at least one parameter mismatch."
      !end if
    else if (restarted .ne. "N") then
      ! The boolean at the end determines whether it reads in a previous forcing file or not
      !call gaussian(numGPcolumns, window_pts, window_skip, timeStop, &
        !tstep, gaussian_tmscl, tol, myid, numTrashLines, .true.)
      call gaussian(t,.true.,.false.)
      print *, "Restoring forcing from dump"
    else
      !call gaussian(numGPcolumns, window_pts, window_skip, timeStop, &
        !tstep, gaussian_tmscl, tol, myid, numTrashLines, .false.)
      call gaussian(t,.false.,.true.)
      print *, "Regenerating a new forcing with initial value zero"
    end if
    ! ------------------------------------------------------------------------------------- 

    
    ! Read Values from the File
    ! ------------------------------------------------------------------------------------- 
    !open(21, file="forcing_data/forcing"//trim(str(myid))//".dat")
      !do i=1, numTrashLines
        !read(21, *)
      !end do
      !do i=1, numGProws
        !read(21,dataFormat) gpTimeVals(i), gpForcingVals(i, :)
        !if(myid .eq. 0) then
        !end if
      !end do
    !close(21)
    ! ------------------------------------------------------------------------------------- 
   
    ! Compute Slope Values for Interpolation 
    ! ------------------------------------------------------------------------------------- 
    !do i=1, numGProws-1
      !interpSlopeVals(i, :) = (gpForcingVals(i+1,:) - gpForcingVals(i,:))/&
        !&(gpTimeVals(i+1)-gpTimeVals(i))
    !end do
    ! ------------------------------------------------------------------------------------- 
  end if
  ! ----------------------------------------------------------------------------------------- 
end subroutine read_stochastic_forcing
