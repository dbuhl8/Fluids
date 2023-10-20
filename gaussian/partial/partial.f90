
program partial

    use gaussian_mod
    
    implicit none

    integer, parameter :: n=20, n2 = 2
    integer, parameter :: dp=selected_real_kind(15)
    real(dp), parameter :: xstep=0.1, sc = 1.0_dp
    integer, parameter :: iter = 100/xstep
    real(dp) :: x(n), y(n), junk, e1, wsc
    real(dp) :: mu=0.0_dp
    real(dp) :: nxy1(2),  xlog(iter)
    character(len = 11) :: tab
    integer :: i, j

    wsc = xstep * n

    
    tab = char(11)
    open(16, file="../plotData/pgen.dat")
    open(13, file="../plotData/opartial.dat")
    
    ! Seeding initial data
    do i = 1, n
        x(i) = (i-1) * xstep
        xlog(i) = x(i)
        y(i) = 0_dp
        write (16, *) xlog(i), tab, y(i)
        write (13, *) xlog(i), tab, y(i)
    end do
        !Calls a subroutine which operates on a window and only updates the current vector after 5-10 iterations
    do i = 1, iter - n
        xlog(i+n) = xlog(i+n-1) + xstep 
        if(mod(i,n2) .eq. 0) then
            call sgnp(x, y, n, xstep * n2, sc)
            write (16, *) xlog(i+n), tab, y(n)
            write (13, *) xlog(i+n), tab, y(n)
        else
            j = mod(i, n2)
            nxy1 = fgnp(x, y, n, xstep * j, sc)
            write (16, *) xlog(i+n), tab, nxy1(2)
        end if
    end do
    close(13)
    close(16)

    !Reading Data to determine where it is centered on
    open(10, file="../plotData/pgen.dat", status="old")
        do i = 1, iter
            read(10,*) junk, tab, e1
            mu = mu + e1
        end do
        mu = mu / iter
    close(10)
    
    print *, " "
    print *, " "
    print *, "OUTPUT LOG FOR PARTIAL.F90"
    print *, " "
    print *, "Number of Timesteps a Gaussian Timescale     : ", xstep/sc
    print *, "Timestep Length                              : ", xstep
    print *, "Y1 centered on                               : ", mu
    print *, " "
    print *, "        ------ Window Details ------         "
    print *, "Window Scale                                 : ", wsc
    print *, "Number of Points in the Window               : ", n
    print *, "The window is updated every                  : ", n2
    print *, "Number of Tao in the Window Length           : ", wsc/sc
    print *, " "
    print *, "The line below is the name of the plot png in the Trails folder (Fluids/plots/trials)"
    print "('Win',i2,'Sc',g3.,'Dx',g3.,'Skip',i1,'.png')", n, sc, xstep, n2
  

end program partial






