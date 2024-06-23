! File: gaussian.f90 formerly partial.f90 from the gaussian directory
! Author: Dante Buhl
! Date: Dec 9, 2023

subroutine gaussian(numGPcols, n, n2, xstop, xstep, sc, tol, id, trashlines, restarting)

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
    use defprecision_module
    
    implicit none

    integer(kind = ki), intent(in) :: n, n2, numGPcols, id ! GP parameters
    real(kind = kr), intent(in) :: xstop, xstep, sc, tol ! GP parameters
    logical, intent(in) :: restarting
    
    integer(kind = ki) :: iter, trashlines, n3, looplen !Iterations
    real(kind = kr) :: x(n) ! Time window
    real(kind = kr) :: y(n, numGPcols) ! Forcing window
    real(kind = kr) :: wsc ! Window Scale
    real(kind = kr), allocatable :: outputs(:,:) !(n2, numGPcols+1) size
    integer :: i, j
    
    iter = NINT(xstop/xstep) !Rounds the dividend between xstop and xstep to an integer
    
    if (numGPcolumns .ne. 0) then

        wsc = xstep * n * n2
        trashlines = 16
        
        open(16, file="forcing_data/forcing"//trim(str(id))//".dat")
        write (16, headerFormatInt)  "Number of Columns generated                  : ", numGPcolumns
        write (16, headerFormatInt)  "Number of Points in the Window               : ", n
        write (16, headerFormatInt)  "The window is updated every                  : ", n2
        write (16, headerFormatReal)  "Gaussian Timescale                           : ", sc
        write (16, headerFormatReal)  "Timestep Length                              : ", xstep
        write (16, headerFormatReal)  "Tolerance                                    : ", tol
        write (16, headerFormatReal)  "Stop at time                                 : ", xstop
        write (16, headerFormatInt)  "Number of lines in header                    : ", trashlines
        write (16, "(A)")  "#"
        write (16, "(A)")  "#       ------ Window Details ------         "
        write (16, headerFormatReal)  "Window Scale                                 : ", wsc
        write (16, headerFormatReal)  "Number of timesteps in the window            : ", wsc / xstep
        write (16, headerFormatReal)  "Window Delta X                               : ", xstep * n2
        write (16, headerFormatReal)  "Number of Tao in the Window Length           : ", wsc/sc
        write (16, headerFormatReal)  "Number of Timesteps in a Gaussian Timescale  : ", sc/xstep
        write (16, "(A)")  "#-------------------------------------------------"


        ! Seeding initial data into the window!
        if(restarting) then
        ! read from file
            open(69, file="forcing_data/restartforcing"//trim(str(id))//".dat")
                print *, "accessing restart forcing file"
                do i = 1, n
                    read(69, dataformat) x(i), y(i, :)
                    write(16, dataformat) x(i), y(i, :)
                end do
            close(69)
        else
            do i = 1, n
                x(i) = (i-1) * xstep
                y(i, :) = 0_kr
                write (16, dataFormat) x(i), y(i, :)
            end do
        end if

        ! This can be adjusted to generate all points between window updates at once
        n3 = mod(iter, n2) !n3 is the number of remaining points that doesn't complete a window (will always be less than n2)
        looplen = (iter-n3)/n2 ! number of times to fill a batch of points and update the window
        allocate(outputs(n2, numGPcols))
        do i = 1, looplen
            outputs = fgnp(x, y, numGPcols, n, n2, xstep, sc, tol) ! generates a batch of new points
            do j = 1, n2
                write (16, dataFormat) outputs(j,:)
            end do              

            ! makes room for a new point, data is translated to the left
            do j = 1, n-1
                x(j) = x(j+1)
                y(j, :) = y(j+1, :)
            end do
           
            ! takes the last point generated and adds it to the data to generate the next batch of points. 
            x(n) = outputs(n2, 1)
            y(n, :) = outputs(n2, 2:numGPcols) 
        end do
        deallocate(outputs)

        ! the last batch of points (not enough to update the window)
        allocate(outputs(n3, numGPcols))
        outputs = fgnp(x, y, numGPcols, n, n3, xstep, sc, tol)
        do j = 1, n3
            write (16, dataFormat) outputs(j,:)
        end do              
        deallocate(outputs)
       
        close(16)
        
    end if

end subroutine gaussian






