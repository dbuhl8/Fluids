! File: nerw_gaussian.f90
! Author: Dante Buhl
! Dependencies: gaussian_mod, forcing_module, defprecision_module
! Date: Dec 9, 2023

subroutine gaussian(t, restarting)
  ! numGPcols is the number of GP needed for that processor
  ! xstop shuold be number of timesteps * .0005
  ! n is the number of timesteps in the gaussian window at any given time. 
  ! n2 is the number of points generated before a new point is added to the window. 
  ! xstop is the last point in the domain that is generated. I.E. if xstop is 100 the GP will go from 0 to 100. 
  ! xstep is the delta x used between GP points. 
  ! sc is the Gaussian Timescale used in the kernal function
  ! tol is the relataive tolerance used in the eigenvalue inversion. Higher tolerance, fewer eigenvalues used. 
  ! id is the processor id that calls this subroutine. Two matching output files are tied to this. 
  ! restarting is a boolean value telling the subroutine to look for previous forcing to restart from

  use gaussian_mod
  use parameter_module
  use defprecision_module
  
  implicit none

  ! new vars
  real(kind=kr) :: next_dump_timestep
  real(kind=kr) :: last_forcing_dump_timestep
  real(kind=kr) :: t
  integer(kind=ki) :: num_new, num_old
  integer(kind=ki) :: num_rows, num_cols
  integer(kind=ki) :: current_time_index
  real(kind=kr),allocatable :: new_points(:), new_x(:)

  ! old vars
  logical, intent(in) :: restarting
 
  ! util vars 
  integer :: i, j
  
  if (FORCING_CPU) then
    ! Seeding initial data into the window!
    if(restarting) then
      ! read from file
      open(101, file="forcing_data/"//fdump_in_file//"_forcing"//&
        trim(str(id))//".dat")
        print *, "accessing restart forcing file"
        read(101, headerFormatInt) num_cols
        read(101, headerFormatInt) num_rows
        ! if statement to make sure num_cols and num_rows match
        if (num_rows .ne. numGProws) then
          print *, "Shape mismatch in array from fdump: number of rows"
        elseif (num_cols .ne. numGPcolumns) then
          print *, "Shape mismatch in array from fdump: number of columns"
        end if
        read(101, headerFormatInt) fidum 
        do i = 1, num_rows
          read(101, dataformat) gpTimevals(i), gpForcingvals(i, :)
        end do
      close(101)
      ! this bit finds the last point in gpTimeVals less than t
      current_time_index = findloc(floor(gpTimevals-t), 0)-1
      restarting = .false.
    elseif(t .eq. 0) then
      gpTimevals(num_rows) = 0.0
      do i = 1, numGProws-1
        gpTimevals(numGProws-i) = -i * tstep
      end do
      gpForcingvals = 0.0
      current_time_index = numGProws-1
      from_scratch = .false.
    else 
      current_time_index = findloc(floor(gpTimevals-t), 0)-1
      restarting = .false.
    end if

    ! need to generate a new batch of points 
    num_new = numGProws-current_time_index
    allocate(new_points(num_new), new_x(num_new))
    new_x(1) = gpTimevals(numGProws) + tstep
    do i = 1, num_new-1
      new_x(i+1) = new_x(i) + tstep
    end do 
    do i = 1, num_cols
      new_points = new_dgr(gpTimeVals, gpForcingvals(:,i), new_x, numGProws,&
         num_new, gauss_scale, tol, idum) 
      ! store new points in forcing columns
      gpForcingvals(1:num_rows-current_time_index+1, i) = &
        gpForcingvals(current_time_index:num_rows, i)
      gpForcingvals(current_time_index+1:,i) = new_points(:)
    end do 
    ! shift old points back, append new points
    gpTimevals(1:numGProws-current_time_index+1) = &
      gpTimevals(current_time_index:numGProws)
    gpTimevals(numGProws-current_time_index+2:) = new_points(:)
  end if
end subroutine gaussian






