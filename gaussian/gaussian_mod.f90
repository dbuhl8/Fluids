!-----------------------------------!
!                                   !
!   File    :   gaussian_mod.f90    !
!                                   !
!   Author  :   Dante Buhl          !
!                                   !
!   Date    :   Oct. 10, 2023       !
!                                   !
!-----------------------------------!

module gaussian_mod


!Static Paramters
real(8), parameter :: pi = atan(1.)*4

!Subroutines and Function definitions
contains 


function dgr(xtrain, ftrain, xsample, n, n2, kscale, char1, char2) result(fsample)

    !----------------------------------------------------------------------------------------------------------------------------!
    ! Parameters                                                                                                                 !
    !----------------------------------------------------------------------------------------------------------------------------!
    ! XTRAIN(:) :   REAL(8)     Vector of x values to train gaussian model on                                                    !
    !                                                                                                                            !
    ! FTRAIN(:) :   REAL(8)     Vector of f values to train gaussian model on                                                    !
    !                                                                                                                            !
    ! XSAMPLE(:):   REAL(8)     Vector of x values to find outputs through the gaussian process                                  !
    !                                                                                                                            !
    ! N, N2     :   INTEGER     As parameters n and n2 are the length of xtrain and xsample respectively.                        !
    !                                                                                                                            !
    ! KSCALE    :   REAL(8)     Length Scale parameter given to the kernel function                                              !
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
    real(8), intent(in) :: xtrain(:), ftrain(:), xsample(:), kscale
    integer, intent(in) :: n, n2
    character, intent(in) :: char1, char2

    !Variables
    integer i, j
    real(8) :: fsample(n2)
    real(8) :: inv_sigma(n, n), sigmak(n, n2), sigmakt(n2, n), sigmakk(n2, n2)
    real(8), allocatable :: Q(:, :), D(:, :)
    
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
        call syed(inv_sigma, n, Q, D)
        Q = transpose(Q)
        fsample = matmul(transpose(matmul(Q, sigmak)), matmul(D, matmul(Q, ftrain)))
        sigmakk = sigmakk - matmul(transpose(matmul(Q, sigmak)), matmul(D, matmul(Q, sigmak)))
    else
        inv_sigma = isqla(inv_sigma, n)
        fsample = matmul(sigmakt, matmul(inv_sigma, ftrain))
        sigmakk = sigmakk - matmul(sigmakt, matmul(inv_sigma, sigmak))
    end if
    sigmakk = sysq(sigmakk, n2) 
    fsample = fsample + matmul(sigmakk, rnorm(n2))

end function dgr



function isqla(A, k) result (B)

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

    integer, intent(in) :: k
    integer :: i, j, lwork1, lwork2, IPIV(k), info=0, kb1, kb2
    real(8), intent(in) :: A(:, :)
    real(8) :: B(k, k)
    real(8), allocatable :: WORK1(:), WORK2(:)
    integer, external :: ilaenv
    
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

end function isqla


function isqed(A, k) result(B)
    
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
    integer, intent(in) :: k
    real(8) :: A(:, :)
    real(8) :: B(k, k), eigen(k)
    real(8), allocatable :: work(:)
    integer :: lwork, kb, info, i
    integer, external :: ilaenv
    
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
    
end function isqed


subroutine syed(A, k, Q, D)

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
    real(8) :: tol = 10**(-2)
    integer, intent(in) :: k
    real(8) :: A(:, :)
    real(8), allocatable :: Q(:, :), D(:, :)
    real(8) :: eigen(k)
    real(8), allocatable :: work(:)
    integer :: lwork, kb, info, i, idx
    integer, external :: ilaenv
    
    kb = ilaenv(1, 'dsytrd', 'L', k, k, -1, -1)
    lwork = (kb + 2) * k

    allocate(work(lwork))

    call dsyev('V', 'L', k, A, k, eigen, work, lwork, info)

    do i = 1, k
        if(eigen(i) .ge. tol) then
            idx = i
            exit
        end if
    end do
    
    allocate(Q(k, (k-idx+1)), D(k-idx+1, k-idx+1))
    D = 0.0
    do i = 1, k-idx+1
        Q(:, i) = A(:, i+idx-1)
        D(i, i) = 1/eigen(i+idx-1)
    end do

end subroutine syed


function sysq(sigma, k)
    
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
    
    real(8), intent(in) :: sigma(:, :)
    integer, intent(in) :: k

    real(8) :: eigen(k), diag(k, k), sysq(k, k)
    integer :: lwork, info, kb, i, j
    real(8), allocatable :: work(:)
    integer, external :: ilaenv

    kb = ilaenv(1, 'dsytrd', 'L', k, k, -1, -1)
    lwork = (kb + 2) * k

    allocate(work(lwork))

    call dsyev('V', 'L', k, sigma, k, eigen, work, lwork, info)

    diag = 0.0
    do i = 1, k
        diag(i, i) = sqrt(eigen(i))
    end do
    sysq = matmul(sigma, diag)

    do j = 1, k
        do i = 1, k
            if (sysq(i, j) .ne. sysq(i, j)) then
                sysq(i, j) = 0
            end if
        end do
    end do

end function sysq


function ek(A, B, na, nb, kscale)
    
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

    real(8), intent(in) :: a(:), b(:), kscale
    integer, intent(in) :: na, nb

    integer :: i, j
    real(8) :: ek(na, nb)
    
    do j = 1,nb
        do i = 1,na
             ek(i, j) = exp(-abs(a(i) - b(j))/kscale)
        end do
    end do

end function ek



function esqk(a, b, na, nb, kscale)
    
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

    real(8), intent(in) :: a(:), b(:), kscale
    integer, intent(in) :: na, nb

    integer :: i, j
    real(8) :: esqk(na, nb)
    
    do j = 1,nb
        do i = 1,na
             esqk(i, j) = exp((-1/(2 * kscale**2)) * (a(i) - b(j))**2)
        end do
    end do

end function esqk


!This function was taken from the Box-Muller Transform wiki page
function rnorm(len)

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
    integer, intent(in) :: len
    real(8) :: vec(len), vec2(len), rnorm(len)
    
    call random_number(vec)
    call random_number(vec2)
    
    rnorm = sqrt(-2 * log(vec)) * cos(2 * pi * vec2)

end function rnorm

function fgnp(x, f, k, delta, sc) result(nxf)

    !----------------------------------------------------------------------------------------------------------------------------!
    ! Parameters                                                                                                                 !
    !----------------------------------------------------------------------------------------------------------------------------!
    ! X(:)      :   REAL(8)     Vector of x coordinates in ascending order, returned with a new point                            !
    !                                                                                                                            !
    ! F(:)      :   REAL(8)     Vector of f coordinates in ascending order, returned with a new point                            !
    !                                                                                                                            !
    ! K         :   REAL(8)     Length of X and F vectors                                                                        !
    !----------------------------------------------------------------------------------------------------------------------------!
    implicit none
    
    integer, intent(in) :: k
    real(8) :: x(:), f(:), sc, delta
    real(8) :: nx(1), nf(1), nxf(2)
    !print *, "x:",x,"f:", f
    nx = x(k) + delta
    nf = dgr(x, f, nx, k, 1, sc, 'S', 'E')
    !print *, "NX:", nx, "NF:", nf 
    nxf(1) = nx(1)
    nxf(2) = nf(1)

end function fgnp

subroutine sgnp(x, f, k, delta, sc)

    !----------------------------------------------------------------------------------------------------------------------------!
    ! Parameters                                                                                                                 !
    !----------------------------------------------------------------------------------------------------------------------------!
    ! X(:)      :   REAL(8)     Vector of x coordinates in ascending order, returned with a new point                            !
    !                                                                                                                            !
    ! F(:)      :   REAL(8)     Vector of f coordinates in ascending order, returned with a new point                            !
    !                                                                                                                            !
    ! K         :   REAL(8)     Length of X and F vectors                                                                        !
    !----------------------------------------------------------------------------------------------------------------------------!
    
    implicit none
    
    integer, intent(in) :: k
    real(8) :: x(:), f(:), sc, delta
    real(8) :: nx(1), nf(1)
    integer :: i

    nx(1) = x(k) + delta
    nf = dgr(x, f, nx, k, 1, sc, 'S', 'E')

    do i = 1, k-1
        x(i) = x(i+1)
        f(i) = f(i+1)
    end do

    x(k) = nx(1)
    f(k) = nf(1)

end subroutine

function np(x, k) result(x2)

    !----------------------------------------------------------------------------------------------------------------------------!
    ! Parameters                                                                                                                 !
    !----------------------------------------------------------------------------------------------------------------------------!
    ! X(:)      :   REAL(8)     Vector of x coordinates in ascending order                                                       !
    !----------------------------------------------------------------------------------------------------------------------------!
    ! X2        :   REAL(8)     New x coordinate                                                                                 !
    !----------------------------------------------------------------------------------------------------------------------------!
    
    implicit none

    real(8), intent(in) :: x(:)
    integer, intent(in) :: k
    
    real(8) :: avgdelta, x2
    integer :: i
    
    do i = 1, k-1
        avgdelta = avgdelta + (x(i+1) - x(i))
    end do
    avgdelta = avgdelta / (k-1)
    x2 = x(k) + avgdelta

end function np


end module gaussian_mod



