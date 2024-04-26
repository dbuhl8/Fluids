
program gen

    use gaussian_mod
    
    implicit none
    integer, parameter :: n=10, init=3, iter=1000
    integer, parameter :: dp=selected_real_kind(15)
    real(8), parameter :: xstep=0.05, sc = 1.25_dp

    real(8) :: x(n), y(n), junk, e1, wsc, nxy(2)
    real(8) :: mu=0.0_dp
    character(len = 11) :: tab
    integer :: i, j

    wsc = xstep * n
    
    print *, "Timescale Non-Dimensionalization Analysis"
    print *, " "
    print *, "Ratio (Timestep / Gaussian Timescale)        : ", xstep/sc
    print *, "Ratio (Timestep / Window Scale)              : ", xstep/wsc
    print *, "Ratio (Gaussian Timescale / Window Timescale): ", sc/wsc


    tab = char(11)
    open(16, file="../plotData/ogen.dat")
    
    ! Seeding initial data
    do i = 1, init
        x(i) = (i-1) * xstep
        y(i) = 0
        !write (16, *) x(i), tab, y(i)
    end do
    
    do i = init, n-1
        nxy = fgnp(x(1:i), y(1:i), i, xstep, sc)
        x(i+1) = nxy(1)
        y(i+1) = nxy(2)
        mu = mu + nxy(2)
        !write (16, *) nxy(1), tab, nxy(2)
    end do
    
    !Calls a subroutine which operates on a window and updates it after each call
    do i = n+1, iter
        if(mod(i,n) .eq. 1) then
            do j = 1, n
                write (16, *) x(j), tab, y(j)
            end do
        end if
        call sgnp(x, y, n, xstep, sc)
        mu = mu + y(n)
    end do

    close(16)
    print *, " "
    print *, "Gaussian Process is centered on              : ", mu/iter
   end program gen






