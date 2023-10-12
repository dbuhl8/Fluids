
program partial

    use gaussian_mod
    
    implicit none

    real(8) :: x(10), y1(10), xstep=0.25, sc = 1.8, junk, e1, e2, e3, e4, kval
    real(8) :: y2(10), y3(10), y4(10), mu1=0., mu2=0., mu3=0., mu4=0.
    real(8) :: x2(10), x3(10), x4(10), nxy1(2), nxy2(2), nxy3(2), nxy4(2), xlog(100)
    real(8), allocatable :: randY(:)
    character(len = 11) :: tab
    integer :: i, j, n=10, iter=100, n2=5

    allocate(randY(n))
    randy = rnorm(n)

    call random_number(kval)

    tab = char(11)
    open(16, file="ogen.dat")
    open(13, file="partial.dat")
    
    ! Seeding initial data
    do i = 1, n
        x(i) = (i-1) * xstep
        xlog(i) = x(i)
        y1(i) = exp(-x(i))
        y2(i) = sin(x(i) * kval)
        y3(i) = cos(x(i) * kval)
        y4(i) = abs(x(i)) - 1
        write (16, *) xlog(i), tab, y1(i), tab, y2(i), tab, y3(i), tab, y4(i)
        write (13, *) xlog(i), tab, y1(i), tab, y2(i), tab, y3(i), tab, y4(i)
    end do
    x2 = x
    x3 = x
    x4 = x
    
    !Calls a subroutine which operates on a window and only updates the current vector after 5-10 iterations
    do i = 1, iter - n
        xlog(i+n) = xlog(i+n-1) + xstep 
        if(mod(i,n2) .eq. 0) then
            do j = 1, n
                call sgnp(x, y1, n, xstep * n2, sc)
                call sgnp(x2, y2, n, xstep * n2, sc)
                call sgnp(x3, y3, n, xstep * n2, sc)
                call sgnp(x4, y4, n, xstep * n2, sc)
                write (16, *) xlog(i+n), tab, y1(n), tab, y2(n), tab, y3(n), tab, y4(n)
            end do
        else
            j = mod(i, n2)
            nxy1 = fgnp(x, y1, n, xstep * j, sc)
            nxy2 = fgnp(x2, y2, n, xstep * j, sc)
            nxy3 = fgnp(x3, y3, n, xstep * j, sc)
            nxy4 = fgnp(x4, y4, n, xstep * j, sc)
            write (16, *) xlog(i+n), tab, nxy1(2), tab, nxy2(2), tab, nxy3(2), tab, nxy4(2)
        end if
        if(mod(i, n2*n) .eq. 0) then
            do j = 1, n
                write (13, *) x(j), tab, y1(j), tab, y2(j), tab, y3(j), tab, y4(j)
            end do
        end if
    end do
    close(13)
    close(16)

    !Reading Data to determine where it is centered on
    open(10, file="ogen.dat", status="old")
        do i = 1, iter
            read(10,*) junk, tab, e1, tab, e2, tab, e3, tab, e4
            mu1 = mu1 + e1
            mu2 = mu2 + e2
            mu3 = mu3 + e3
            mu4 = mu4 + e4
        end do
        mu1 = mu1 / iter
        mu2 = mu2 / iter
        mu3 = mu3 / iter
        mu4 = mu4 / iter
        print *, "Y1 centered on :", mu1
        print *, "Y2 centered on :", mu2
        print *, "Y3 centered on :", mu3
        print *, "Y4 centered on :", mu4
        print *, " "
        print *, "Ratio (Timestep / Gaussian Timescale) : ", xstep/sc
        print *, "Ratio (Timestep / Period)             : ", xstep*kval/(2*pi)
        print *, "Ratio (Gaussian Timescale / Period)   : ", sc*kval/(2*pi)
    close(10)

end program partial






