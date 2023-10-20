
program partial

    use gaussian_mod
    
    implicit none

    integer, parameter :: n=20, n2 = 2
    integer, parameter :: dp=selected_real_kind(15)
    real(dp), parameter :: xstep=0.1, sc = 1.0_dp
    integer, parameter :: iter = 100/xstep
    real(dp) :: x(n), y1(n), junk, e1, wsc
    real(dp) :: mu1=0.0_dp
    real(dp) :: nxy1(2),  xlog(iter)
    character(len = 11) :: tab
    integer :: i, j

    wsc = xstep * n

    print *, "Timescale Non-Dimensionalization Analysis"
    print *, " "
    print *, "Ratio (Timestep / Gaussian Timescale)        : ", xstep/sc
    print *, "Ratio (Timestep / Window Scale)              : ", xstep/wsc\
    print *, "------ Window Details ------"
    print *, "Number of points in the window at a given time" 
    print *, "Number of Tao (Gaussian Timescales) in the Window Length: ", wsc/sc
  

    tab = char(11)
    open(16, file="../plotData/pgen.dat")
    open(13, file="../plotData/opartial.dat")
    
    ! Seeding initial data
    do i = 1, n
        x(i) = (i-1) * xstep
        xlog(i) = x(i)
        y1(i) = 0
        write (16, *) xlog(i), tab, y1(i)
        write (13, *) xlog(i), tab, y1(i)
    end do
        !Calls a subroutine which operates on a window and only updates the current vector after 5-10 iterations
    do i = 1, iter - n
        xlog(i+n) = xlog(i+n-1) + xstep 
        if(mod(i,n2) .eq. 0) then
            call sgnp(x, y1, n, xstep * n2, sc)
            write (16, *) xlog(i+n), tab, y1(n)
            write (13, *) xlog(i+n), tab, y1(n)
        else
            j = mod(i, n2)
            nxy1 = fgnp(x, y1, n, xstep * j, sc)
            write (16, *) xlog(i+n), tab, nxy1(2)
        end if
    end do
    close(13)
    close(16)

    !Reading Data to determine where it is centered on
    open(10, file="../plotData/pgen.dat", status="old")
        do i = 1, iter
            read(10,*) junk, tab, e1
            mu1 = mu1 + e1
        end do
        mu1 = mu1 / iter
        print *, "Y1 centered on :", mu1
        
    close(10)

end program partial






