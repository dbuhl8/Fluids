!File: shearSolve.f90
!Author: Dante Buhl
!Dependencies: LinAl.mod, LAPACK, LBLAS

program shearSolve

  use LinAl

  implicit none
  include 'mpif.h'

  integer, parameter :: M = 10, N = 10 !"Num of Fourier Modes in X (N) and Y(M)"
  integer, parameter :: LDA = 5*(2*N+1)*(2*M+1) !"Leading Dimension of A"
  real :: ky = 1, kx = 0.5 !"Wavenumbers for the Fourier Decomp"
  real :: kzmin = 0.0000, kzmax = 10 !"Range of KZ values"

  ! Indices for the run
  integer :: i, j, k, indu, indv, indw, indt, indp, indm, l

  ! Non-Dimensional quantities
  real :: Pe, Re, Ri, f, dispterm

  ! Coefficients for the fourier decomp
  real :: v0 = 0.915,  alpha=1.0
  complex, parameter :: cmplx_1=(1.0, 0.0), cmplx_i=(0.0, 1.0)

  ! LAPACK Routine Vars
  character :: jobvl = "N", jobvr = "V"
  complex, dimension(LDA) :: ALFA, BETA
  complex, allocatable :: A(:, :), B(:, :), V(:, :), VL(:, :)
  complex, allocatable :: VR(:, :)
  complex, allocatable :: work(:)
  real, dimension(8*LDA) :: RWORK=0.0
  integer :: lwork, info, nk
  real, allocatable :: D(:, :)
  ! This is delcared as real because we will return only the real part of the
  ! eigenvalue

  ! Variables for measuring compute time
  real :: start, finish

  ! Variables for collective output
  real, allocatable :: veclambda(:, :), veckz(:, :)
  real :: kz, lambda

  ! MPI variables
  integer :: myid, ie, np
  integer :: msgid, src, dest
  integer :: buffer
  integer :: stat(MPI_STATUS_SIZE)
  character(9) :: signal


  ! starting MPI
  call MPI_INIT(ie)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ie)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)
  call MPI_BARRIER(MPI_COMM_WORLD,ie)

  ! set number of wavenumbers to number of processors
  nk = np
  
  ! This determines which values for kz, this processor will run on
  kz = kzmin + (myid)*(kzmax - kzmin)/(np-1)
  print "(A, I3, A, F8.4)", "Processor: ", myid, ", Kz: ", kz
  ! declaring the non-dimensional parameters for this run
  Pe = 10000.0
  Re = 10000.0
  Ri = 100.0
  ! The floquet coefficient is programmed to work but is not actually used
  f = 0.0

  ! Allocating memory for the LAPACK routine
  allocate(A(LDA, LDA), B(LDA, LDA), D(LDA, LDA), V(LDA, LDA), VL(LDA, LDA), &
          VR(LDA, LDA))

  A = 0.0
  B = 0.0
  ! this is used to see how long each eigensolve takes
  call cpu_time(start)

  do k = -M, M
    do i = -N, N
      ! computing the indices for this loop
      indu = (i + N + 1) + (2*N+1)*(M+k)
      indv = (2*N + 1)*(2*M + 1) + indu
      indw = (2*N + 1)*(2*M + 1) + indv
      indt = (2*N + 1)*(2*M + 1) + indw
      indp = (2*N + 1)*(2*M + 1) + indt

      dispterm = ((f+i)**2)*kx**2 + ((f+k)**2)*ky**2 + kz**2
      ! check if i > -N
      if (i > -N) then ! q_(n-1, m)
        !---------------------------------------------------------
        !w equation
        A(indw, indw-1) = (f+k)*ky*alpha*v0*0.5*cmplx_i !
        !---------------------------------------------------------
        !t equation
        A(indt, indt-1) = (f+k)*ky*alpha*v0*0.5*cmplx_i !
        !---------------------------------------------------------
        !v equation
        A(indv, indu-1) = kx*alpha*v0*0.5*cmplx_i  !
        A(indv, indv-1) = (f+k)*ky*alpha*v0*0.5*cmplx_i  !
        !---------------------------------------------------------
        !u equation
        A(indu, indu-1) = (f+k)*ky*alpha*v0*0.5*cmplx_i 
        !---------------------------------------------------------

        ! check if k > -M (i > -N)
        if (k > -M) then !q_(n-1, m-1)
          !---------------------------------------------------------------
          !w equation
          A(indw, indw-2*N-2) = alpha*0.5*(-kx*(f+i-1)+ky*(f+k-1)*0.5)*cmplx_i !
                                !-(i+f-1)*kx*alpha*0.5*(0.0, 1.0) - &
                                !(f+k-1)*ky*alpha*v1*0.5*(0.0, 1.0)
          !---------------------------------------------------------------
          !t equation
          A(indt, indt-2*N-2) = alpha*0.5*(-kx*(f+i-1)+ky*(f+k-1)*0.5)*cmplx_i !
                                !-(i+f-1)*kx*alpha*0.5*(0.0, 1.0) - &
                                !(f+k-1)*ky*alpha*v1*0.5*(0.0, 1.0) 
          !---------------------------------------------------------------
          !v equation
          A(indv, indv-2*N-2) = alpha*0.5*(-kx(f+i-1)+ky*(f+k)*0.5)*cmplx_i !
                                !-(i+f-1)*kx*alpha*0.5*(0.0, 1.0) -&
                                !ky*alpha*v1*0.5*(0.0, 1.0) - &
                                !((f+k-1)*ky*alpha*v1*0.5)*(0.0, 1.0)
          A(indv, indu-2*N-2) = (kx*alpha*0.25)*cmplx_i !
          !---------------------------------------------------------------
          !u equation
          A(indu, indu-2*N-2) = alpha*0.5*(-(f+i)*kx+ky*(f+k+1)*0.5)*cmplx_i !
                                !-(kx*alpha*0.5)*(0.0, 1.0) - &
                                 !(f+k-1)*ky*alpha*v1*0.5*(0.0, 1.0) -&
                                 !(f+i-1)*kx*alpha*0.5*(0.0, 1.0)
          A(indu, indv-2*N-2) = -ky*alpha*0.5*cmplx_i !
          !---------------------------------------------------------------
        end if
        ! check if k < M (i > -N)
        if (k < M) then !q_(n-1, m+1)
          !-------------------------------------------------------------------
          ! w equation
          A(indw, indw+2*N) = alpha*0.5*(-kx*(f+i-1)-ky*(f+k+1)*0.5)*cmplx_i !
                              !-(i+f-1)*kx*alpha*0.5*(0.0, 1.0) + &
                              !(f+k+1)*ky*alpha*v_1*0.5*(0.0, 1.0) 
          !-------------------------------------------------------------------
          ! t equation
          A(indt, indt+2*N) = alpha*0.5*(-kx*(f+i-1)-ky*(f+k+1)*0.5)*cmplx_i !
                              !-(i+f-1)*kx*alpha*0.5*(0.0, 1.0) + &
                               !(f+k+1)*ky*alpha*v_1*0.5*(0.0, 1.0) 

          !-------------------------------------------------------------------
          ! v equation
          A(indv, indv+2*N) = alpha*0.5*(-kx*(f+i-1)-(f+k)*ky*0.5)*cmplx_i !
                              !-(i+f-1)*kx*alpha*0.5*(0.0, 1.0) + &
                              !(f+k+1)*ky*alpha*v_1*0.5*(0.0, 1.0) +&
                               !(ky*alpha*v_1*0.5)*(0.0, 1.0)
          A(indv, indu+2*N) =  -kx*alpha*0.25*cmplx_i !

          !-------------------------------------------------------------------
          ! u equation
          A(indu, indu+2*N) = alpha*0.5*(-kx*(f+i)-ky*(f+k+1)*0.5)*cmplx_i !
                              !-(kx*alpha*0.5)*(0.0, 1.0) + &
                               !(f+k+1)*ky*alpha*v_1*0.5*(0.0, 1.0) -&
                               !(f+i-1)*kx*alpha*0.5*(0.0, 1.0)
          A(indu, indv+2*N) =  ky*alpha*0.5*cmplx_i !
          !-------------------------------------------------------------------
        end if
      end if

      ! this logic section is not needed yet for this problem
      ! check if i < N
      if(i < N) then ! q_(n+1, m)
        !---------------------------------------------------------
        ! w equation
        A(indw, indw+1) = (f+k)*ky*alpha*v0*0.5*cmplx_i !
        !---------------------------------------------------------
        ! t equation
        A(indt, indt+1) = (f+k)*ky*alpha*v0*0.5*cmplx_i !
        !---------------------------------------------------------
        ! u equation
        A(indu, indu+1) = (f+k)*ky*alpha*v0*0.5*cmplx_i !
        !---------------------------------------------------------
        ! v equation
        A(indv, indv+1) = (f+k)*ky*alpha*v0*0.5*cmplx_i !
        A(indv, indu+1) = -kx*alpha*v0*0.5*cmplx_i !
        !---------------------------------------------------------
        ! check if k > -M (i < N)
        if (k > -M) then ! q_(n+1, m-1) 
          !-------------------------------------------------------------------
          ! w equation
          A(indw, indw-2*N) = alpha*0.5*(-kx*(f+i+1)-ky*(f+k-1)*0.5)*cmplx_i !
                              !-(i+f+1)*kx*alpha*0.5*(0.0, 1.0) + &
                               !(f+k-1)*ky*alpha*v_1*0.5*(0.0, 1.0) 
          !-------------------------------------------------------------------
          ! t equation
          A(indt, indt-2*N) = alpha*0.5*(-kx*(f+i+1)-ky*(f+k-1)*0.5)*cmplx_i !
                              !-(i+f+1)*kx*alpha*0.5*(0.0, 1.0) +&
                               !(f+k-1)*ky*alpha*v_1*0.5*(0.0, 1.0) 
          !-------------------------------------------------------------------
          ! u equation
          A(indu, indu-2*N) = alpha*0.5*(-(f+i)*kx+ky*(f+k-1)*0.5)*cmplx_i !
                              !kx*alpha*0.5*(0.0, 1.0) - &
                             !(f+i+1)*kx*alpha*0.5*(0.0, 1.0) &
                              !+ (f+k-1)*ky*alpha*v1*0.5*(0.0, 1.0)
          A(indu, indv-2*N) = -ky*alpha*0.5*cmplx_i !
          !-------------------------------------------------------------------
          ! v equation
          A(indv, indv-2*N) = alpha*0.5*(-kx*(f+i+1)-ky*(f+k)*0.5)*cmplx_i !
                               !-ky*alpha*v_1*0.5*(0.0, 1.0) - &
                               !(f+i+1)*kx*alpha*0.5*(0.0, 1.0) &
                               !- (f+k-1)*ky*alpha*v_1*0.5*(0.0, 1.0)
          A(indv, indu-2*N) = kx*alpha*0.25*cmplx_i !
          !-------------------------------------------------------------------
        end if
        ! check if k < M (i < N)
        if (k < M) then ! q_(n+1, m+1)
          !--------------------------------------------------------------------
          ! w equation
          A(indw,indw+2*N+2) = alpha*0.5*(-kx*(f+i+1)+ky*(f+k+1)*0.5)*cmplx_i !
                                 !-(i+f+1)*kx*alpha*0.5*(0.0, 1.0) -&
                                 !(f+k+1)*ky*alpha*v_1*0.5*(0.0, 1.0) 
          !--------------------------------------------------------------------
          ! t equation
          A(indt,indt+2*N+2) = alpha*0.5*(-kx*(f+i+1)+ky*(f+k+1)*0.5)*cmplx_i !
                               !-(i+f+1)*kx*alpha*0.5*(0.0, 1.0) -&
                                 !(f+k+1)*ky*alpha*v_1*0.5*(0.0, 1.0) 
          !--------------------------------------------------------------------
          ! u equation
          A(indu, indu+2*N+2) = alpha*0.5*(-kx*(f+i)+ky*(f+k+1)*0.5)*cmplx_i !
                                !kx*alpha*0.5*(0.0, 1.0) - & 
                                !(f+i+1)*kx*alpha*0.5*(0.0, 1.0) & 
                                !- (f+k+1)*ky*alpha*v1*0.5*(0.0, 1.0)
          A(indu, indv+2*N+2) = ky*alpha*0.5*cmplx_i !
          !--------------------------------------------------------------------
          ! v equation
          A(indv, indv+2*N+2) = alpha*0.5*(-kx*(f+i+1)+ky*(f+k)*0.5)*cmplx_i !
                               !-ky*alpha*v_1*0.5*(0.0, 1.0) - &
                                !(f+i+1)*kx*alpha*0.5*(0.0, 1.0) & 
                              !+ (f+k+1)*ky*alpha*v_1*0.5*(0.0, 1.0)
          A(indv, indu+2*N+2) = -kx*alpha*0.25*cmplx_i !
          !--------------------------------------------------------------------
        end if
      end if 

      if (k > -M) then ! q_(n, m-1)
        !---------------------------------------------------------
        ! w equation
        A(indw, indw-(2*N+1)) = -(f+i)*kx*0.5*cmplx_1 !
        !---------------------------------------------------------
        ! t equation
        A(indt, indt-(2*N+1)) = -(f+i)*kx*0.5*cmplx_1 !
        !---------------------------------------------------------
        ! v equation
        A(indv, indv-(2*N+1)) = -(f+i)*kx*0.5*cmplx_1 !
        !---------------------------------------------------------
        ! u equation
        A(indu, indu-(2*N+1)) = -(f+i)*kx*0.5*cmplx_1 !
        A(indu, indv-(2*N+1)) = -ky*0.5*cmplx_1 !
        !---------------------------------------------------------
      end if
      ! check if k < M
      if (k < M) then ! q_(n, m+1)
        !---------------------------------------------------------
        ! w equation
        A(indw, indw+(2*N+1)) = (f+i)*kx*0.5*cmplx_1 !
        !---------------------------------------------------------
        ! t equation
        A(indt, indt+(2*N+1)) = (f+i)*kx*0.5*cmplx_1 !
        !------------------------------------------------------
        ! v equation
        A(indv, indv+(2*N+1)) = (f+i)*kx*0.5*cmplx_1 !
        !------------------------------------------------------
        ! u equation
        A(indu, indu+(2*N+1)) = (f+i)*kx*0.5*cmplx_1 !
        A(indu, indv+(2*N+1)) = -ky*0.5*cmplx_1 !
        !---------------------------------------------------------
      end if
     
      ! Dissipitive Terms
      A(indu,indu)= -dispterm/Re
      A(indw,indw)= -dispterm/Re
      A(indv,indv)= -dispterm/Re
      A(indt,indt)= -dispterm/Pe
      
      ! Stratification Terms here
      A(indt,indw)=-1.0
      A(indw,indt)= Ri

      ! Pressure Terms
      A(indu,indp)=-kx*(f+i)  ! avoid complex matrix by defining variable ip
      A(indv,indp)=-ky*(f+k)
      A(indw,indp)=-kz

      ! Eigenvalue Matrix B, 1 on diagonal
      B(indu,indu)=1.0
      B(indv,indv)=1.0
      B(indw,indw)=1.0
      B(indt,indt)=1.0
      B(indp,indp)=0.0
      
      ! Continuity Equation
      A(indp,indu) = kx*(f+i)
      A(indp,indv) = ky*(f+k)
      A(indp,indw) = kz
    end do
  end do

  info = 0
  VL = 0.0
  VR = 0.0
  ALFA = 0.0
  BETA = 0.0
  lwork = -1
  allocate(work(1))
  work = (0.0, 0.0)
  ! Arguments are numbered to help with debugging, (LAPACK will return error
  ! messages with the argument number that is incorrect
  !           1       2     3   4   5   6   7    8     9    10   11  12  13  
  call cggev(jobvl, jobvr, LDA, A, LDA, B, LDA, ALFA, BETA, VL, LDA, VR, LDA,&
                                                     WORK, LWORK, RWORK, info)
  !                                                   14     15    16     17
  print "(A, I3, A, 2F10.2)", "Processor: ",myid,", Computed lwork: ",work(1)
  lwork = int(work(1))
  deallocate(work)  
  allocate(work(lwork))
  !           1       2     3   4   5   6   7    8     9    10   11  12  13 
  call cggev(jobvl, jobvr, LDA, A, LDA, B, LDA, ALFA, BETA, VL, LDA, VR, LDA,&
                                                     WORK, LWORK, RWORK, info)
  !                                                   14     15    16     17
  deallocate(work)

  if (info .ne. 0) then
    print *, "Error in LAPACK Routine"
  end if
  do l = 1, LDA
    if(abs(BETA(l)) > 10.d-6) then
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
  print "(A, I3, A, F8.3)", "Processor: ", myid, ", Computed Lambda: ", lambda
  print "(A, I3, A, F8.3, A)", "Processor: ", myid, ", Done with eigenproblem&
                                   . Time elapsed: ", finish-start, " seconds"

  !bring all of the lambda's together call MPI_GATHER(lambda, numeigen,
  !MPI_REAL8, veclambda, nk, MPI_REAL8, 0, MPI_COMM_WORLD) call MPI_GATHER(kz,
  !numeigen, MPI_REAL8, veckz, nk, MPI_REAL8, 0, MPI_COMM_WORLD)

  !if on cpu 0
  !if (myid .eq. 0) then
      !print *, "Eigenvalues"
      !call printmat(veclambda, nk, 1)
      !print *, "Kz values"
      !call printmat(veckz, nk, 1)

      !open(10, file="kz.dat")
          !do i = 1, nk
          !write (10, *) veckz(i, 1) 
          !end do
      !close(10)
      !open(12, file="lambda.dat")
          !do i = 1, nk
              !write (12, *) veclambda(i, 1)
          !end do
      !close(12)
      !open(13, file="nk.dat")
          !write (13, *) nk
      !close(13)
      !deallocate(veclambda, veckz)
  !end if

  deallocate(A, B, V, VL, D, VR)

  call MPI_FINALIZE(ie)

end program shearSolve




