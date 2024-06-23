!File: Gaussian Mod
!Author: Dante Buhl
!Dependencies: defprecision_module, parameter_module, other module, MPI,


module gaussian_mod
  use defprecision_module
  use parameter_module, only : pi
  use state_module, only : decomp_independent_random

  implicit none
    
  ! public vars
  integer :: idx,idy,idz

  !Subroutines and Function definitions
  contains 
    function dgr(xtrain, ftrain, xsample, n, n2, kscale, tol, char1, char2) result(fsample)
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Parameters                                                                                                                 !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! XTRAIN(:) :   REAL(8)     Vector of x values to train gaussian model on                                                    !
      !                                                                                                                            !
      ! FTRAIN(:) :   REAL(8)     Vector of f values to train gaussian model on                                                    !
      !
      ! NUMCOLS   :   INTEGER     Num of columns in FTRAIN
      !                                                                                                                            !
      ! XSAMPLE(:):   REAL(8)     Vector of x values to find outputs through the gaussian process                                  !
      !                                                                                                                            !
      ! N, N2     :   INTEGER     As parameters n and n2 are the length of xtrain and xsample respectively.                        !
      !                                                                                                                            !
      ! KSCALE    :   REAL(8)     Length Scale parameter given to the kernel function                                              !
      !                                                                                                                            !
      ! TOL       :   REAL(8)     Relative Tolerance for eigenvalue retention in Eigendecomposition
      !                                                                                                                            !
      ! Char1     :   CHAR        Determines whether Expontial or Exponential Squared Kernel is used. Set to 'S' for Exp. Squared, !
      !                           anything else for Exp.                                                                           !
      !                                                                                                                            !
      ! Char2     :   CHAR        Determins whether eigendecomposition or block diagonal factorization in order to compute inverse !
      !                           matrices. Set to 'E' for Eigendecomposition, anything else for block diagonal factorization      !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Output                                                                                                                     !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! FSAMPLE(:):   REAL(8)     Vector of f values found through the gaussian process                                            !
      !----------------------------------------------------------------------------------------------------------------------------!
      implicit none
      
      !Intent IN
      real(kind=kr), intent(in) :: xtrain(:), ftrain(:), xsample(:), kscale, tol
      integer(kind = ki), intent(in) :: n, n2
      character, intent(in) :: char1, char2

      !Variables
      integer(kind=ki) :: i, j, k
      real(kind=kr) :: fsample(n2)
      real(kind=kr) :: inv_sigma(n, n), sigmak(n, n2), sigmakt(n2, n), sigmakk(n2, n2)
      real(kind=kr), allocatable :: Q(:, :), D(:, :)
      
      if(char1 .eq. 'S') then     
        inv_sigma = esqk(xtrain, xtrain, n, n, kscale)
        sigmak = esqk(xtrain, xsample, n, n2, kscale)
        sigmakt = transpose(sigmak)
        sigmakk = esqk(xsample, xsample, n2, n2, kscale)
      else 
        inv_sigma = ek(xtrain, xtrain, n, n, kscale)
        sigmak = ek(xtrain, xsample, n, n2, kscale)
        sigmakt = transpose(sigmak)
        sigmakk = ek(xsample, xsample, n2, n2, kscale)
      end if
      if(char2 .eq. 'E') then
        call syed(inv_sigma, n, Q, D, tol)
        Q = transpose(Q)
        fsample = matmul(transpose(matmul(Q, sigmak)), matmul(D, matmul(Q, ftrain)))
        sigmakk = sigmakk - matmul(transpose(matmul(Q, sigmak)), matmul(D, matmul(Q, sigmak)))
      else
        inv_sigma = isqla(inv_sigma, n)
        fsample = matmul(sigmakt, matmul(inv_sigma, ftrain))
        sigmakk = sigmakk - matmul(sigmakt, matmul(inv_sigma, sigmak))
      end if
      sigmakk = sysq(sigmakk, n2)
      ! mu_star = mu + sigma_star * Normal distribution
      fsample = fsample + matmul(sigmakk, rnorm(n2))
      deallocate(Q, D)
    end function dgr

    function new_dgr(old_x,old_f,new_x,num_old,num_new,gauss_scale,tol,idum)
      ! Produces a gaussian process regression based off of (old_x, old_f)
      implicit none
      ! IO Vars
      real(kind=kr) :: old_x(:), old_f(:), new_x(:), gauss_scale, tol
      integer(kind=ki) :: num_old, num_new, idum
      real(kind=kr), dimension(num_new) :: new_dgr

      ! Lin Al Vars
      real(kind=kr), dimension(num_old,num_old) :: Koo
      real(kind=kr), dimension(num_old,num_new) :: Kon
      real(kind=kr), dimension(num_new,num_old) :: Kno
      real(kind=kr), dimension(num_new,num_new) :: Knn
      real(kind=kr), allocatable :: Q(:,:), D(:,:)

      Koo = esqk(old_x, old_x, num_old, num_old, gauss_scale)
      Kon = esqk(old_x, new_x, num_old, num_new, gauss_scale)
      Kno = transpose(Kon)
      Knn = esqk(new_x, new_x, num_new, num_new, gauss_scale)

      call syed(koo, num_old, Q, D, tol)
      Q = transpose(Q)
      new_dgr = matmul(transpose(matmul(Q, Kon)), matmul(D, matmul(Q, old_f)))
      Knn = Knn-matmul(transpose(matmul(Q, Kon)),matmul(D,matmul(Q,Kon)))
      
      knn = sysq(knn, num_new)
      ! mu_star = mu + sigma_star * Normal distribution
      new_dgr = new_dgr + matmul(knn, rnorm(num_new,idum))
      deallocate(Q, D)
    end function new_dgr

    !generates numRows amount of new points for a gaussian process with training set x, f
    function fgnp(x, f, numCols, k, numRows, delta, sc, tol) result(newxfvals)
      !            (x, f, numCols, k, numRows, delta, sc, tol)
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Parameters                                                                                                                 !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! X(:)      :   REAL(8)     Vector of x coordinates in ascending order, returned with a new point                            !
      !                                                                                                                            !
      ! F(:, :)      :   REAL(8)     Vector of f coordinates in ascending order, returned with a new point                            !
      !                                                                                                                            !
      ! K         :   REAL(8)     Length of X and F vectors                                                                        !
      !----------------------------------------------------------------------------------------------------------------------------!
      implicit none
     
      integer(kind=ki), intent(in) :: numCols , numRows
      integer(kind=ki), intent(in) :: k
      integer(kind=ki) :: i
      real(kind=kr) :: x(:), f(:, :), sc, delta, tol
      real(kind=kr) :: newxvals(numRows), newfvals(numRows), newxfvals(numRows, numcols+1)
    
      do i = 1, numRows 
        newxvals(i) = x(k) + delta*i
      end do
      newxfvals(:, 1) = newxvals(:)
      
      do i = 1, numCols
        newfvals = dgr(x, f(:, i), newxvals, k, numRows, sc, tol, 'S','E')
        newxfvals(:, i+1) = newfvals(:)
      end do
    end function fgnp

    subroutine sgnp(x, f, numCols, k, delta, sc, tol)
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Parameters                                                                                                                 !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! X(:)      :   REAL(8)     Vector of x coordinates in ascending order, returned with a new point                            !
      !                                                                                                                            !
      ! F(:, :)      :   REAL(8)     Vector of f coordinates in ascending order, returned with a new point                            !
      !                                                                                                                            !
      ! K         :   REAL(8)     Length of X and F vectors                                                                        !
      !----------------------------------------------------------------------------------------------------------------------------!
      
      implicit none
     
      integer(kind=ki), intent(in) :: numCols 
      integer(kind=ki), intent(in) :: k
      real(kind=kr) :: x(:), f(:, :), sc, delta, tol
      real(kind=kr) :: npx(1), npf(1), npfholder(numCols)
      integer(kind=ki) :: i, j

      npx(1) = x(k) + delta   
      do j = 1, numCols
        npf = dgr(x, f(:, j), npx, k, 1, sc, tol, 'S', 'E')
        npfholder(j) = npf(1)
      end do

      do i = 1, k-1
        x(i) = x(i+1)
        f(i, :) = f(i+1, :)
      end do

      x(k) = npx(1)
      f(k, :) = npfholder(:)
    end subroutine sgnp

    function isqla(A, k) result (B)
      ! ISQLA - Inverse of a SQuare symmetric matrix LApack version
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Parameters                                                                                                                 !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! A(:,:)    :   REAL(8)     Square Symmetric Matrix to be inverted                                                           !
      !                                                                                                                            !
      ! K         :   INTEGER     Length of A on both dimensions                                                                   !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Output                                                                                                                     !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! B(:)      :   REAL(8)     The computed inverse of A                                                                        !
      !----------------------------------------------------------------------------------------------------------------------------!
      
      implicit none

      integer(kind=ki), intent(in) :: k
      integer(kind=ki) :: i, j, lwork1, lwork2, IPIV(k), info=0, kb1, kb2
      real(kind = kr), intent(in) :: A(:, :)
      real(kind = kr) :: B(k, k)
      real(kind = kr), allocatable :: WORK1(:), WORK2(:)
      integer(kind=ki), external :: ilaenv
      
      kb1 = ilaenv(1, 'dsytrf', 'L', k, k, -1, -1)
      kb2 = ilaenv(1, 'dsytri2', 'L', k, k, -1, -1)
      
      lwork1 = k * kb1
      lwork2 = (k + kb2 + 1) * (kb2 + 3)

      allocate(WORK1(lwork1), WORK2(lwork2))

      B = A   

      call dsytrf('L', k, B, k, IPIV, WORK1, lwork1, info)
      call dsytri2('L', k, B, k, IPIV, WORK2, lwork2, info)
      do j = 1, k
        do i = 1, k
          if(j > i) then
            B(i, j) = B(j, i)
          end if
        end do
      end do

      deallocate(WORK1, WORK2)
    end function isqla

    function isqed(A, k) result(B)
      !ISQED - Inverse of a SQuare symmetric matrix using EigenDecomposition
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Parameters                                                                                                                 !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! A(:,:)    :   REAL(8)     Square Symmetric Matrix to be inverted                                                           !
      !                                                                                                                            !
      ! K         :   INTEGER     Length of A on both dimensions                                                                   !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Output                                                                                                                     !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! B(:)      :   REAL(8)     The computed inverse of A                                                                        !
      !----------------------------------------------------------------------------------------------------------------------------!
      
      implicit none
      real(8) :: tol = 10**(-2)
      integer(kind=ki), intent(in) :: k
      real(kind = kr) :: A(:, :)
      real(kind = kr) :: B(k, k), eigen(k)
      real(kind = kr), allocatable :: work(:)
      integer(kind=ki) :: lwork, kb, info, i
      integer(kind=ki), external :: ilaenv
      
      kb = ilaenv(1, 'dsytrd', 'L', k, k, -1, -1)
      lwork = (kb + 2) * k

      allocate(work(lwork))

      call dsyev('V', 'L', k, A, k, eigen, work, lwork, info)

      tol = tol * eigen(k)
      do i = 1, k
        if(eigen(i) .le. tol) then
          B(i, i) = 0.0
          A(:, i) = 0.0
        else 
          B(i, i) = 1/eigen(i)
        end if
      end do
      B = matmul(A, matmul(B, transpose(A)))

      deallocate(work)
    end function isqed

    subroutine syed(A, k, Q, D, tol)
      !SYED - SYmetric EigenDecomposition
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Parameters                                                                                                                 !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! A(:,:)    :   REAL(8)     Square Symmetric Matrix to be inverted                                                           !
      !                                                                                                                            !
      ! K         :   INTEGER     Length of A on both dimensions                                                                   !
      !                                                                                                                            !
      ! Q(:,:)    :   REAL(8)     Unallocated Double array, will be allocated in subroutine                                        !
      !                                                                                                                            !
      ! D(:,:)    :   REAL(8)     Unallocated Double array, will be allocated in subroutine
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Output                                                                                                                     !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Q(:,:)    :   REAL(8)     Matrix of eigenvectors which correspond to the eigenvalues in D                                  !
      !                                                                                                                            !
      ! D(:,:)    :   REAL(8)     Inv Diagonal Matrix of Eigenvalues such that no eigenvalues are below the tolerance              !
      !----------------------------------------------------------------------------------------------------------------------------!
      
      implicit none
      real(kind = kr), intent(in) :: tol
      integer(kind=ki), intent(in) :: k
      real(kind = kr) :: A(:, :)
      real(kind = kr), allocatable :: Q(:, :), D(:, :)
      real(kind = kr) :: eigen(k)
      real(kind = kr), allocatable :: work(:)
      integer(kind=ki) :: lwork, kb, info, i, j
      !integer(kind=ki), external :: ilaenv

      ! I can just call the function twice      
      !kb = ilaenv(1, 'dsytrd', 'L', k, k, -1, -1)
      !lwork = (kb + 2) * k
      allocate(work(1))
      lwork = -1
      call dsyev('V', 'U', k, A, k, eigen, work, lwork, info)
      lwork = work(1)
      deallocate(work)
      allocate(work(lwork))
      call dsyev('V', 'U', k, A, k, eigen, work, lwork, info)
      deallocate(work)

      if (info .ne. 0) then
        print *, "Error in LAPACK Routine in SYED routine"
      end if

      ! This routine returns an array of eigenvalues in ascending order

      ! this constructs the psuedo inverse matrices Q, and D
      do i = 1, k
        if(eigen(i) .ge. tol*abs(eigen(k))) then
          j = i
          exit
        end if
      end do

      allocate(Q(k,(k-j+1)), D(k-j+1,k-j+1))
      D = 0.0
      do i = j, k
        Q(:, i-j+1) = A(:, i)
        D(i-j+1,i-j+1) = 1./eigen(i)
      end do
    end subroutine syed

    function sysq(sigma, k)
      !SYSQ - SYmmetric SQuare root of a matrix 
      !---------------------------------------------------------------!
      !   Parameters                                                  !
      !---------------------------------------------------------------!
      ! SIGMA(:,:)    :   REAL(8) Positive Defintite Symmetric Matrix !
      !                                                               !
      ! K             :   INTEGER Number of Rows/Columns in SIGMA     !
      !---------------------------------------------------------------!
      !   Output                                                      !
      !---------------------------------------------------------------!
      ! SYSQ(:,:)     :   REAL(8) "Square-Root" of SIGMA              !
      !---------------------------------------------------------------!

      implicit none
      
      real(kind = kr), intent(in) :: sigma(:, :)
      integer(kind=ki), intent(in) :: k

      real(kind = kr) :: eigen(k), diag(k, k), sysq(k, k)
      integer(kind=ki) :: lwork, info, kb, i, j
      real(kind = kr), allocatable :: work(:)
      integer, external :: ilaenv

      !kb = ilaenv(1, 'dsytrd', 'L', k, k, -1, -1)
      !lwork = (kb + 2) * k
      allocate(work(1))
      lwork = -1
      call dsyev('V', 'U', k, sigma, k, eigen, work, lwork, info)
      lwork = work(1)
      deallocate(work)
      allocate(work(lwork))
      call dsyev('V', 'U', k, sigma, k, eigen, work, lwork, info)
      deallocate(work)

      if (info .ne. 0) then
        print *, "Error in LAPACK Routine in SYSQ routine"
      end if

      diag = 0.0
      do i = 1, k
        diag(i, i) = sqrt(eigen(i))
      end do
      sysq = matmul(sigma, diag)

      ! NOTE: Commenting this out to see if somethign bad happens 
      ! why is this here?
      !do j = 1, k
        !do i = 1, k
          !if (sysq(i, j) .ne. sysq(i, j)) then
            !sysq(i, j) = 0
          !end if
        !end do
      !end do
      ! it zeros the lower triangular matrix where it is not symmetric?
      deallocate(work)
    end function sysq

    function ek(A, B, na, nb, kscale)
      !EK - Exponential Kernel
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Parameters                                                                                                                 !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! A(:)      :   REAL(8)     Vector of doubles                                                                                !
      !                                                                                                                            !
      ! B(:)      :   REAL(8)     Vector of doubles                                                                                !
      !                                                                                                                            !
      ! NA/NB     :   INTEGER     Length of A/B vectors respectively                                                               !
      !                                                                                                                            !
      ! KSCALE    :   REAL(8)     Scaling parameter of the kernel                                                                  !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Output                                                                                                                     !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! EK(:,:)   :   REAL(8)     Exponential Kernel Matrix of DIM(NA, NB)                                                         !
      !----------------------------------------------------------------------------------------------------------------------------!

      implicit none

      real(kind = kr), intent(in) :: a(:), b(:), kscale
      integer(kind=ki), intent(in) :: na, nb

      integer(kind=ki) :: i, j
      real(kind = kr) :: ek(na, nb)
      
      do j = 1,nb
        do i = 1,na
          ek(i, j) = exp(-abs(a(i) - b(j))/kscale)
        end do
      end do
    end function ek

    function esqk(a, b, na, nb, kscale)
      !ESQK - Exponential SQuared Kernel 
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Parameters                                                                                                                 !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! A(:)      :   REAL(8)     Vector of doubles                                                                                !
      !                                                                                                                            !
      ! B(:)      :   REAL(8)     Vector of doubles                                                                                !
      !                                                                                                                            !
      ! NA/NB     :   INTEGER     Length of A/B vectors respectively                                                               !
      !                                                                                                                            !
      ! KSCALE    :   REAL(8)     Scaling parameter of the kernel                                                                  !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Output                                                                                                                     !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! ESQK(:,:) :   REAL(8)     Exponential Kernel Matrix of DIM(NA, NB)                                                         !
      !----------------------------------------------------------------------------------------------------------------------------!

      implicit none

      real(kind = kr), intent(in) :: a(:), b(:), kscale
      integer(kind=ki), intent(in) :: na, nb

      integer(kind=ki) :: i, j
      real(kind = kr) :: esqk(na, nb)
      
      do j = 1,nb
        do i = 1,na
          esqk(i, j) = exp((-1/(2 * kscale**2)) * (a(i) - b(j))**2)
        end do
      end do
    end function esqk

    !This function was taken from the Box-Muller Transform wiki page
    ! This needs to be modified to be deterministically random
    function rnorm(len,idum)
      !RNORM - Random NORMal distribution sample (len number of independent samples)
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Parameters                                                                                                                 !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! LEN       :   INTEGER     Length of returned Normal Distribution Sample                                                    !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! Output                                                                                                                     !
      !----------------------------------------------------------------------------------------------------------------------------!
      ! RNORM(:)  :   REAL(8)     Vector of samples from Normal Distribution                                                       !
      !----------------------------------------------------------------------------------------------------------------------------!

      implicit none
      integer(kind=ki), intent(in) :: len, idum
      real(kind = kr) :: vec(len), vec2(len), rnorm(len), rn

      ! this still needs some work to give it a unique seed
      do i = 1, len
        vec(i) = ran2(idum) !forcing_decomp_independent_random(params, idum)
        vec2(i) = ran2(idum) !forcing_decomp_independent_random(params, idum)
      end do 
      
      rnorm = sqrt(-2. * log(vec)) * cos(2. * pi * vec2)
    end function rnorm

    function forcing_decomp_independent_random(params, idum)
      implicit none
      integer(kind=ki) :: idum
      ! params
      !return ran2(some seed here)
    end function forcing_decomp_independent_random
end module gaussian_mod



