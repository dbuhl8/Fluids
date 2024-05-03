!File: shearSolve.f90
!Author: Dante Buhl
!Dependencies: MPI, LAPACK, LBLAS

program shearSolve

  implicit none
  include 'mpif.h'

  integer, parameter :: kr=kind(dble(0.))
  integer, parameter :: Nmax=5, Mmax=5,nk=10 !"Num of Fourier Modes in X and Y" 
  integer, parameter :: LDA = 5*(2*Nmax+1)*(2*Mmax+1) !"Leading Dimension of A"
  real(kind=kr) :: ky = 1, kx = 0.5 !"Wavenumbers for the Fourier Decomp"
  real(kind=kr) :: kzmin = 10.d-10, kzmax = 10 !"Range of KZ values"
  real(kind=kr) :: Delta_Kz=0.05, dkz=0.05 !Step distasnce between kz
  logical :: bool = .TRUE.

  ! Indices for the run
  integer :: i, j, k, indu, indv, indw, indt, indp, indm, l, n, m

  ! Non-Dimensional quantities
  real(kind=kr) :: Pe, Re, Ri, f, dispterm

  ! Coefficients for the fourier decomp
  real(kind=kr) :: v0 = 0.914539,  alpha=1.0
  complex(kind=kr), parameter :: cmplx_1=(1.0, 0.0), cmplx_i=(0.0, 1.0)

  ! LAPACK Routine Vars
  character :: jobvl = "N", jobvr = "V"
  complex(kind=kr), dimension(LDA) :: ALFA, BETA
  complex(kind=kr), allocatable :: A(:, :), B(:, :), V(:, :), VL(:, :)
  complex(kind=kr), allocatable :: VR(:, :)
  complex(kind=kr), allocatable :: work(:)
  real(kind=kr), dimension(8*LDA) :: RWORK=0.0
  integer :: lwork, info
  real(kind=kr), allocatable :: D(:, :)
  ! This is delcared as real because we will return only the real part of the
  ! eigenvalue

  ! Variables for measuring compute time
  real(kind=kr) :: start, finish

  ! Variables for collective output
  real(kind=kr), dimension(nk) :: lambda_r, kz_r
  real(kind=kr), allocatable :: kz_p(:), lambda_p(:)
  real(kind=kr) :: kz, lambda

  ! MPI variables
  integer :: myid, ie, np
  integer :: msgid, src, dest
  integer :: buffer
  integer, allocatable :: kz_per_proc(:)
  integer, allocatable :: disp_vec(:)
  integer :: stat(MPI_STATUS_SIZE)
  character(9) :: signal
  integer :: randomthing

  ! starting MPI
  call MPI_INIT(ie)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ie) ! get this processor id
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie) ! get number of processors
  call MPI_BARRIER(MPI_COMM_WORLD,ie) ! I have absolutely no idea what this does

  allocate(kz_per_proc(np), disp_vec(np))
  ! Determine how many kz's to run on each proccessor
  ! if nk/np doesn't divide evenly, then some processors will get more
  ! eigenvalue problems to solve than others. This gode will load the extra
  ! solves on the last processors. 
  kz_per_proc = nk/np
  if (mod(nk, np) .ne. 0) then
    do i = 0, mod(nk, np)-1
      kz_per_proc(np-i) = kz_per_proc(np-i) + 1
    end do 
  end if

  lambda_r = 0.0
  kz_r = 0.0

  allocate(kz_p(kz_per_proc(myid+1)), lambda_p(kz_per_proc(myid+1)))
  kz_p = 0.0
  lambda_p = 0.0
  if (myid .eq. 0) then
    ! Calculating delta_kz 
    if(bool) then 
      Delta_Kz = (kzmax - kzmin)/(nk-1)
    else 
      Delta_kz = dkz !here dkz is a default value for Delta_Kz
    end if
    do i = 1, nk
      kz_r(i) = kzmin + (i-1)*Delta_kz
    end do
  end if

  ! Distribute kz values to each processor
  disp_vec = 0
  do i = 1, np-1
    disp_vec(i+1) = sum(kz_per_proc(1:i))
  end do
                     !sbuf,  scounts,   displs,    stype,  rbuf, 
  call MPI_SCATTERV(kz_r, kz_per_proc, disp_vec, MPI_REAL8, kz_p,&
                    kz_per_proc(myid+1), MPI_REAL8, 0, MPI_COMM_WORLD, ie)
                     !rcount,             rtype,   root    comm       ierr
  lambda_p = 0.0

  ! Declaring the non-dimensional parameters for this run
  ! this can be dynamically programmed with a fortran namelist 
  Pe = 10000.0
  Re = 10000.0
  Ri = 100.0
  ! The floquet coefficient is programmed to work but is not actually used
  f = 0.0

  ! Allocating memory for the LAPACK routine
  allocate(A(LDA, LDA), B(LDA, LDA), D(LDA, LDA), V(LDA, LDA), VL(LDA, LDA), &
          VR(LDA, LDA))


  ! DO PER KZ, this solves the eigenvalue problem
  do i = 1, kz_per_proc(myid+1)
    kz = kz_p(i)
    A = 0.0
    B = 0.0

    ! this is used to see how long each eigensolve takes
    call cpu_time(start)

    do m = -Mmax, Mmax
      do n = -Nmax, Nmax
        ! computing the indices for this loop
        indu = (n + Nmax + 1) + (2*Nmax+1)*(Mmax+m)
        indv = (2*Nmax + 1)*(2*Mmax + 1) + indu
        indw = (2*Nmax + 1)*(2*Mmax + 1) + indv
        indt = (2*Nmax + 1)*(2*Mmax + 1) + indw
        indp = (2*Nmax + 1)*(2*Mmax + 1) + indt

        !Note: one exclamation mark indicates its been checked against my algebra,
        !and two also checked against pascale's code.

        ! check if i > -Nmax
        if (n > -Nmax) then ! q_(n-1, m)
          !---------------------------------------------------------
          !w equation
          A(indw, indw-1) = (f+m)*ky*alpha*v0*0.5*cmplx_i !!
          !---------------------------------------------------------
          !t equation
          A(indt, indt-1) = (f+m)*ky*alpha*v0*0.5*cmplx_i !!
          !---------------------------------------------------------
          !v equation
          A(indv, indu-1) = kx*alpha*v0*0.5*cmplx_i  !!
          A(indv, indv-1) = (f+m)*ky*alpha*v0*0.5*cmplx_i  !!
          !---------------------------------------------------------
          !u equation
          A(indu, indu-1) = (f+m)*ky*alpha*v0*0.5*cmplx_i !!
          !---------------------------------------------------------

          ! check if k > -M (i > -Nmax)
          if (m > -Mmax) then !q_(n-1, m-1)
            !---------------------------------------------------------------
            !w equation
            A(indw,indw-2*Nmax-2) = alpha*0.5*(-kx*(f+n-1)+ky*(f+m-1)*0.5)*cmplx_i!!
            !---------------------------------------------------------------
            !t equation
            A(indt,indt-2*Nmax-2) = alpha*0.5*(-kx*(f+n-1)+ky*(f+m-1)*0.5)*cmplx_i!!
            !---------------------------------------------------------------
            !v equation
            A(indv,indv-2*Nmax-2) = alpha*0.5*(-kx*(f+n-1)+ky*(f+m)*0.5)*cmplx_i!!
            A(indv,indu-2*Nmax-2) = (kx*alpha*0.25)*cmplx_i !!
            !---------------------------------------------------------------
            !u equation
            A(indu,indu-2*Nmax-2) = alpha*0.5*(-(f+n)*kx+ky*(f+m-1)*0.5)*cmplx_i!!
            A(indu,indv-2*Nmax-2) = -ky*alpha*0.5*cmplx_i !!
            !---------------------------------------------------------------
          end if
          ! check if m < M (i > -Nmax)
          if (m < Mmax) then !q_(n-1, m+1)
            !-------------------------------------------------------------------
            ! w equation
            A(indw, indw+2*Nmax) = alpha*0.5*(-kx*(f+n-1)-ky*(f+m+1)*0.5)*cmplx_i!!
            !-------------------------------------------------------------------
            ! t equation
            A(indt, indt+2*Nmax) = alpha*0.5*(-kx*(f+n-1)-ky*(f+m+1)*0.5)*cmplx_i!!
            !-------------------------------------------------------------------
            ! v equation
            A(indv,indv+2*Nmax) = alpha*0.5*(-kx*(f+n-1)-(f+m)*ky*0.5)*cmplx_i !!
            A(indv,indu+2*Nmax) =  -kx*alpha*0.25*cmplx_i !!
            !-------------------------------------------------------------------
            ! u equation
            A(indu,indu+2*Nmax) = alpha*0.5*(-kx*(f+n)-ky*(f+m+1)*0.5)*cmplx_i !!
            A(indu,indv+2*Nmax) =  ky*alpha*0.5*cmplx_i !!
            !-------------------------------------------------------------------
          end if
        end if

        ! check if n < Nmax
        if(n < Nmax) then ! q_(n+1, m)
          !---------------------------------------------------------
          ! w equation
          A(indw, indw+1) = (f+m)*ky*alpha*v0*0.5*cmplx_i !!
          !---------------------------------------------------------
          ! t equation
          A(indt, indt+1) = (f+m)*ky*alpha*v0*0.5*cmplx_i !!
          !---------------------------------------------------------
          ! u equation
          A(indu, indu+1) = (f+m)*ky*alpha*v0*0.5*cmplx_i !!
          !---------------------------------------------------------
          ! v equation
          A(indv, indv+1) = (f+m)*ky*alpha*v0*0.5*cmplx_i !!
          A(indv, indu+1) = -kx*alpha*v0*0.5*cmplx_i !!
          !---------------------------------------------------------

          ! check if m > -Mmax (n < Nmax)
          if (m > -Mmax) then ! q_(n+1, m-1) 
            !-------------------------------------------------------------------
            ! w equation
            A(indw, indw-2*Nmax) = alpha*0.5*(-kx*(f+n+1)-ky*(f+m-1)*0.5)*cmplx_i!!
            !-------------------------------------------------------------------
            ! t equation
            A(indt, indt-2*Nmax) = alpha*0.5*(-kx*(f+n+1)-ky*(f+m-1)*0.5)*cmplx_i!!
            !-------------------------------------------------------------------
            ! u equation
            A(indu, indu-2*Nmax) = alpha*0.5*(-(f+n)*kx-ky*(f+m-1)*0.5)*cmplx_i !!
            A(indu, indv-2*Nmax) = -ky*alpha*0.5*cmplx_i !!
            !-------------------------------------------------------------------
            ! v equation
            A(indv, indv-2*Nmax) = alpha*0.5*(-kx*(f+n+1)-ky*(f+m)*0.5)*cmplx_i !!
            A(indv, indu-2*Nmax) = kx*alpha*0.25*cmplx_i !!
            !-------------------------------------------------------------------
          end if
          ! check if m < Mmax (n < Nmax)
          if (m < Mmax) then ! q_(n+1, m+1)
            !--------------------------------------------------------------------
            ! w equation
            A(indw,indw+2*Nmax+2) = alpha*0.5*(-kx*(f+n+1)+ky*(f+m+1)*0.5)*cmplx_i!!
            !--------------------------------------------------------------------
            ! t equation
            A(indt,indt+2*Nmax+2) = alpha*0.5*(-kx*(f+n+1)+ky*(f+m+1)*0.5)*cmplx_i!!
            !--------------------------------------------------------------------
            ! u equation
            A(indu,indu+2*Nmax+2) = alpha*0.5*(-kx*(f+n)+ky*(f+m+1)*0.5)*cmplx_i!!
            A(indu,indv+2*Nmax+2) = ky*alpha*0.5*cmplx_i !!
            !--------------------------------------------------------------------
            ! v equation
            A(indv,indv+2*Nmax+2) = alpha*0.5*(-kx*(f+n+1)+ky*(f+m)*0.5)*cmplx_i!!
            A(indv,indu+2*Nmax+2) = -kx*alpha*0.25*cmplx_i !!
            !--------------------------------------------------------------------
          end if
        end if 

        ! check if m > -Mmax
        if (m > -Mmax) then ! q_(n, m-1)
          !---------------------------------------------------------
          ! w equation
          A(indw, indw-(2*Nmax+1)) = -(f+n)*kx*0.5*cmplx_1 !!
          !---------------------------------------------------------
          ! t equation
          A(indt, indt-(2*Nmax+1)) = -(f+n)*kx*0.5*cmplx_1 !!
          !---------------------------------------------------------
          ! v equation
          A(indv, indv-(2*Nmax+1)) = -(f+n)*kx*0.5*cmplx_1 !!
          !---------------------------------------------------------
          ! u equation
          A(indu, indu-(2*Nmax+1)) = -(f+n)*kx*0.5*cmplx_1 !!
          A(indu, indv-(2*Nmax+1)) = -ky*0.5*cmplx_1 !!
          !---------------------------------------------------------
        end if

        ! check if m < Mmax
        if (m < Mmax) then ! q_(n, m+1)
          !---------------------------------------------------------
          ! w equation
          A(indw, indw+(2*Nmax+1)) = (f+n)*kx*0.5*cmplx_1 !!
          !---------------------------------------------------------
          ! t equation
          A(indt, indt+(2*Nmax+1)) = (f+n)*kx*0.5*cmplx_1 !!
          !------------------------------------------------------
          ! v equation
          A(indv, indv+(2*Nmax+1)) = (f+n)*kx*0.5*cmplx_1 !!
          !------------------------------------------------------
          ! u equation
          A(indu, indu+(2*Nmax+1)) = (f+n)*kx*0.5*cmplx_1 !!
          A(indu, indv+(2*Nmax+1)) = -ky*0.5*cmplx_1 !!
          !---------------------------------------------------------
        end if
       
        ! Dissipitive Terms
        dispterm = ((f+n)**2)*kx**2 + ((f+m)**2)*ky**2 + kz**2
        A(indu,indu)= -dispterm/Re !!
        A(indw,indw)= -dispterm/Re !!
        A(indv,indv)= -dispterm/Re !!
        A(indt,indt)= -dispterm/Pe !!
        
        ! Stratification Terms here
        A(indt,indw)=-1.0 !!
        A(indw,indt)= Ri !!

        ! Pressure Terms
        A(indu,indp)=-kx*(f+n) !!
        A(indv,indp)=-ky*(f+m) !!
        A(indw,indp)=-kz !!

        ! Eigenvalue Matrix B, 1 on diagonal
        B(indu,indu)=1.0 !!
        B(indv,indv)=1.0 !!
        B(indw,indw)=1.0 !!
        B(indt,indt)=1.0 !!
        B(indp,indp)=0.0 !!
        
        ! Continuity Equation
        A(indp,indu) = kx*(f+n) !!
        A(indp,indv) = ky*(f+m) !!
        A(indp,indw) = kz !!
      end do
    end do

    ! Initialize to zero for good measure
    info = 0
    VL = 0.0
    VR = 0.0
    ALFA = 0.0
    BETA = 0.0
    lwork = -1
    allocate(work(1))
    work = (0.0, 0.0)
    ! The first call is used to determine the optimal size of the work array
     !           1       2     3   4   5   6   7    8     9    10   11  12  13  
    call zggev(jobvl, jobvr, LDA, A, LDA, B, LDA, ALFA, BETA, VL, LDA, VR, LDA,&
                                                       WORK, LWORK, RWORK, info)
    !                                                   14     15    16     17
    lwork = int(work(1))!Int is needed here since work is a complex-valued array
    deallocate(work)  
    allocate(work(lwork))
    ! Arguments are numbered to help with debugging. LAPACK returns error
    ! messages correlated to the argument number
    !           1       2     3   4   5   6   7    8     9    10   11  12  13 
    call zggev(jobvl, jobvr, LDA, A, LDA, B, LDA, ALFA, BETA, VL, LDA, VR, LDA,&
                                                       WORK, LWORK, RWORK, info)
    !                                                   14     15    16     17
    deallocate(work)

    if (info .ne. 0) then
      print *, "Error in LAPACK Routine"
    end if
    do l = 1, LDA
      if(abs(BETA(l)) > 10.d-14) then
        D(l, l) = real(ALFA(l)/BETA(l))
      end if
    end do
    indm = 1 
    do l = 1, LDA
      if(D(l, l) .le. 0) then
        D(l, l) = 0.;
      end if
      if(D(l, l) .ge. 100.0) then
        D(l, l) = 0.;
      end if
      if(D(l, l) > D(indm, indm)) then
        indm = l
      end if
    end do
    lambda = D(indm, indm)
    call cpu_time(finish)
    
    print "('#',A, I3, 3(A, F12.4), A)", "Processor: ",myid,", Kz :",kz,&
                              ", Computed Lambda :",lambda,  &
                              ", Time elapsed: ",finish-start," seconds"
    lambda_p(i) = lambda
  end do 
  ! END DO PER KZ

  ! Gather computed lambdas into root processor
                     !sbuf,        scounts,            stype,  rbuf
  call MPI_GATHERV(lambda_p, kz_per_proc(myid+1), MPI_REAL8, lambda_r,&
                   kz_per_proc, disp_vec, MPI_REAL8, 0, MPI_COMM_WORLD, ie)
                     !rcounts,   displs,  rtype,   root, comm,        err

  ! Write to OUT file. 
  if (myid .eq. 0) then
    open(15, file='OUT1')
    write(15, "('#', A)")       "      Sinshear results        "
    write(15, "('#', A)")       "------------------------------"
    write(15, "('#', A, I6)")   "Number of Processor's used   :", np
    write(15, "('#', A, I6)")   "Number of Fourier Modes in X :", 2*Nmax+1
    write(15, "('#', A, I6)")   "Number of Fourier Modes in Y :", 2*Mmax+1
    write(15, "('#', A, F10.3)")"Kx                           :", kx
    write(15, "('#', A, F10.3)")"Ky                           :", ky
    write(15, "('#', A, F10.3)")"Reynolds   Number            :", Re
    write(15, "('#', A, F10.3)")"Reynolds   Number            :", Re
    write(15, "('#', A, F10.3)")"Peclet     Number            :", Pe
    write(15, "('#', A, F10.3)")"Richardson Number            :", Ri
    write(15, "('#', A, F10.3)")"Floquet Coefficient          :", f
    write(15, "('#', A, F10.3)")"A                            :", alpha
    do i = 1, nk
      write(15, "(2F20.12)") kz_r(i), lambda_r(i)
    end do 
    close(15)
  end if

  deallocate(A, B, V, VL, D, VR, kz_p, lambda_p, kz_per_proc, disp_vec)

  call MPI_FINALIZE(ie)

  contains 

    subroutine makemode(eigenvec, kz_loc, u)
      ! Note: This subroutine is within the scope of the rest of the program.
      ! Therefore we have that all of the vars declared in shearSolve are
      ! present in makemode. This is why some variables are not passed as
      ! arguments and others are written with "_loc" for a locally scoped
      ! variant. 

      implicit none
      ! Subroutine Arguments
      real(kind=kr) :: eigenvec(:), kz_loc
      ! NOTE: Derived Datatype "velocity" will be copied from PADDI later.
      !type(velocity) :: u
      real(kind=kr) :: u  

      ! Domain Vars
      real(kind=kr), parameter :: pi=cos(-1.0)
      real(kind=kr) :: dx, dy, dz, fourier_mode
      real(kind=kr) :: x, y, z
      integer,parameter :: nx=100, ny=50, nz=50
      
      ! Utility Vars
      integer :: i_loc, j_loc, k_loc, n_loc, m_loc
      integer :: indu_loc, indv_loc, indw_loc

      ! Velocity and Vorticity
      real,dimension(nx, ny, nz) :: ux, uy, uz
      real,dimension(nx, ny, nz) :: wx, wy, wz
      real :: sum_ux, sum_uy, sum_uz
      real :: sum_wx, sum_wy, sum_wz
      real :: unm, vnm, wnm

      ! Initialize variables
      dx = 4.*pi/nx
      dy = 2.*pi/ny
      dz = 2.*pi/nz

      ux = 0.0
      uy = 0.0
      uz = 0.0
      wx = 0.0
      wy = 0.0
      wz = 0.0
      ! Loop through all points in the domain
      do k_loc = 1, nz
        z = k_loc*dz
        do j_loc = 1, ny
          y = j_loc*dy
          do i_loc = 1, nx
            x = i_loc*dx
            sum_ux = 0.0
            sum_uy = 0.0
            sum_uz = 0.0
            sum_wx = 0.0
            sum_wy = 0.0
            sum_wz = 0.0
            ! Compute sum for each component
            do n_loc = -Nmax, Nmax
              do m_loc = -Mmax, Mmax
                ! Find local indices
                indu_loc = (n_loc + Nmax + 1) + (2*Nmax+1)*(Mmax+m_loc)
                indv_loc = (2*Nmax + 1)*(2*Mmax + 1) + indu_loc
                indw_loc = (2*Nmax + 1)*(2*Mmax + 1) + indv_loc
                !NOTE to ask Pascale, why is kz not included here
                fourier_mode = exp(cmplx_i*(n_loc*kx*x+m_loc*ky*y))
                unm =eigenvec(indu_loc)*fourier_mode
                vnm =eigenvec(indv_loc)*fourier_mode
                wnm =eigenvec(indw_loc)*fourier_mode
                sum_ux = sum_ux + unm
                sum_uy = sum_uy + vnm
                sum_uz = sum_uz + wnm
                sum_wx = sum_wx + (cmplx_i*m*k*wnm - &
                             cmplx_i*kz*vnm)*fourier_mode
                sum_wy = sum_wy + (cmplx_i*kz*unm - &
                             cmplx_i*n*kx*wnm)*fourier_mode
                sum_wz = sum_wz + (cmplx_i*n*kx*vnm - &
                              cmplx_i*m*ky*unm)*fourier_mode
              end do 
            end do 
            fourier_mode = exp(cmplx_i*kz*z)
            sum_ux = sum_ux*fourier_mode
            sum_uy = sum_uy*fourier_mode           
            sum_uz = sum_uz*fourier_mode
            sum_wx = sum_wx*fourier_mode
            sum_wy = sum_wy*fourier_mode
            sum_wz = sum_wz*fourier_mode
            ux(i, j, k) = real(sum_ux)
            uy(i, j, k) = real(sum_uy)
            uz(i, j, k) = real(sum_uz)
            wx(i, j, k) = real(sum_wx)
            wy(i, j, k) = real(sum_wy)
            wz(i, j, k) = real(sum_wz)
          end do 
        end do 
      end do
    end subroutine  makemode

end program shearSolve

