! File: gaussian.f90 formerly partial.f90 from the gaussian directory
! Author: Dante Buhl
! Date: Dec 9, 2023

    ! numGPcols is the number of GP needed for that processor
    ! xstop shuold be number of timesteps * .0005
    ! n is the number of timesteps in the gaussian window at any given time. 
    ! n2 is the number of points generated before a new point is added to the window. 
    ! xstop is the last point in the domain that is generated. I.E. if xstop is 100 the GP will go from 0 to 100. 
    ! xstep is the delta x used between GP points. 
    ! sc is the Gaussian Timescale used in the kernal function
    ! tol is the relataive tolerance used in the eigenvalue inversion. Higher tolerance, fewer eigenvalues used. 
    ! id is the processor id that calls this subroutine. Two matching output files are tied to this. 

subroutine gaussian(numGPcols, n, n2, xstop, xstep, sc, tol, id, trashlines)
    use gaussian_mod
    use defprecision_module
    
    implicit none

    integer(kind = ki), intent(in) :: n, n2, numGPcols, id ! GP parameters
    real(kind = kr), intent(in) :: xstop, xstep, sc, tol ! GP parameters
    
    integer(kind = ki) :: iter, trashlines!Iterationsm
    real(kind = kr) :: x(n) ! X window
    real(kind = kr) :: y(n, numGPcols) ! Y window
    real(kind = kr) :: wsc ! Window Scale
    real(kind = kr) :: nxy1(numGPcols+1) ! Point put in dat file but not in the window
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

        print *, "cpu "//trim(str(id))//": got past header write"

        ! Seeding initial data into the window!
        do i = 1, n
            x(i) = (i-1) * xstep
            y(i, :) = 0_kr
            write (16, dataFormat) x(i), y(i, :)
        end do

        print *, "cpu "//trim(str(id))//": got past first point write"

        do i = n, iter
            if(mod(i,n2) .eq. 0) then
                ! Subroutine Gaussian New Point
                call sgnp(x, y, numGPcols, n, xstep * n2, sc, tol) !sgnp updates the window
                write (16, dataFormat) x(n), y(n, :) ! write new point to the dat file
            else
                j = mod(i, n2)
                ! Function Gaussian New Point
                nxy1 = fgnp(x, y, numGPcols, n, xstep * j, sc, tol) !fgnp doesn't update the window
                write (16, dataFormat) nxy1(:) ! write new point to the dat file
            end if
        end do

        close(16)
        
    end if

    print *, "cpu "//trim(str(myid))//": gaussian processes made"

    contains
    
    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

end subroutine gaussian






