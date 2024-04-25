module LinAl 

implicit none 
integer, save :: msize, nsize 
real, dimension(:,:), allocatable, save :: mat
contains

!********************************************************

subroutine readMat(filename)

    implicit none
    character(len=*) :: filename

    integer :: i,j

    open(10,file=filename)

    ! Read the matrix dimensions
    read(10,*) i,j

    ! Read matrix
    do i=1,msize
       read(10,*) ( mat(i,j), j=1,nsize )
    end do

    close(10)
    
  end subroutine readMat

  subroutine trace(A, m, tr)
    
    implicit none
    
    real, intent(in) :: A(:, :)
    integer, intent(in) :: m
    real, intent(out) :: tr 
    integer :: i

    ! sums along the diagonal
    do i = 1, m
        tr = tr + A(i, i)
    end do

  end subroutine trace

  subroutine twonorm(v, nrm)
    
    implicit none

    real, intent(in) :: v(:)
    real, intent(out) :: nrm 

    !square roots the sum of each vector element squared
    nrm = sqrt(sum(v**2))

  end subroutine twonorm

  subroutine printmat(A, m, n)

    real, intent(in) :: A(:, :)
    integer, intent(in) :: m, n
    integer :: i, j

    print "(A, I3, A, I3)", "This is a ", m," by ", n," matrix."
    do i = 1, m
        print "("//trim(str(n))//"F10.3)", A(i, :)
    end do

  end subroutine printmat

  character(len=20) function str(k)
  ! "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str

  subroutine GE(A, B, ma, nb, bool)
        
    implicit none
    
    integer, intent(in) :: ma, nb
    logical :: bool
    real :: A(:, :), B(:, :)
    real, allocatable :: rowapivot(:), rowbpivot(:)
    real, dimension(ma, ma) :: L
    integer :: i, j, k
    bool = .false. 
    
    !initializing our matrices. 
    L = 0.0
    do i = 1, ma
        L(i, i) = 1.0 
    end do

    !begin GE 
    do i = 1, ma-1
        !pivoting step
        k = i
        do j = i+1, ma
            if (abs(A(j, i)) .gt. abs(A(k, i))) then
                k = j
            end if
        end do
        allocate(rowapivot(k:ma), rowbpivot(nb))
        rowapivot = A(i, i:ma)
        rowbpivot = B(i, :)
        A(i, i:ma) = A(k, i:ma)
        B(i, :) = B(k, :)
        A(k, i:ma) = rowapivot
        B(k, :) = rowbpivot
        deallocate(rowapivot, rowbpivot)
        if (A(i, i) .eq. 0) then
            bool = .true.
            return
        end if
        !L Step
        L(i+1:ma, i) = A(i+1:ma,i)/A(i, i)
    
        A(i+1:ma, i:ma) = A(i+1:ma, i:ma) - matmul(L(i+1:ma, i:i), A(i:i, i:ma))
        do j = i+1, ma
            B(j, :) = B(j, :) - L(j, i)*B(i, :)
        end do
    end do 
    
  end subroutine GE

  subroutine backsub(U, B, X, mu, mb)
    implicit none
    
    integer, intent(in) :: mu, mb
    real, intent(in) :: B(:, :), U(:, :)
    real :: X(:, :)
    integer :: i, j, k

    X = 0.0
    X(mu, :) = B(mu, :)/U(mu,mu)
    do i = 1, mu-1
        k = mu - i
        X(k, :) = (B(k, :) - matmul(U(k, k+1:mu), X(k+1:mu, :)))/U(k, k)
    end do

  end subroutine backsub

  subroutine LU(A, ma, bool, P)
    
    implicit none
    
    integer, intent(in) :: ma
    logical :: bool
    real :: A(:, :)
    real, allocatable :: rowapivot(:)
    real, allocatable :: rowlpivot(:)
    real :: rowppivot
    real, dimension(ma, ma) :: L
    integer :: P(:)
    integer :: i, j, k
    bool = .false. 
    !initializing our matrices. 
    L = 0.0
    P = 0
    do i = 1, ma
        P(i) = i
    end do
    k = 0
    !begin GE 
    do i = 1, ma-1
        !pivoting step
        k = i
        do j = i+1, ma
            if (abs(A(j, i)) .gt. abs(A(k, i))) then !set pivot row to k
                k = j
            end if
        end do
        if (k .ne. i) then !pivot A and P
            allocate(rowapivot(k:ma))
            rowapivot = A(i, i:ma)
            rowppivot = P(i)
            A(i, i:ma) = A(k, i:ma)
            P(i) = P(k)
            A(k, i:ma) = rowapivot
            P(k) = rowppivot
            deallocate(rowapivot)
            if (i .ne. 1) then !pivot l
                allocate(rowlpivot(1:i-1))
                rowlpivot = L(i, 1:i-1)
                L(i, 1:i-1) = L(k, 1:i-1)
                L(k, 1:i-1) = rowlpivot
                deallocate(rowlpivot)
            end if
        end if
        if (A(i, i) .eq. 0) then !return singular
            bool = .true.
            return
        end if
        !L U step
        L(i+1:ma, i) = A(i+1:ma,i)/A(i, i)
        do  j = i+1, ma
            A(j, i:ma) = A(j, i:ma) - L(j, i)*A(i, i:ma)
        end do
    end do 
!    print *, "Printing L before adding L to A"
!    call printmat(L, ma, ma)
    A = A + L !this gives A as L + U where the only overlap is along the diagonals

  end subroutine LU

  subroutine LUsolve(LU, ma, B, X, mb, P)

    implicit none

    integer, intent(in) :: ma, mb, P(:)
    real, intent(in) :: LU(:, :), X(:,:)
    real :: B(:, :)
    real, dimension(ma, ma) :: L, U
    real, dimension(ma, mb) :: Y
    integer :: i , j, k

    L = 0.
    U = 0.

    do j = 1, ma
        do i = 1, ma
            if(i .eq. j) then   
                L(i, j) = 1.
                U(i, j) = LU(i, j)
            else if(i .gt. j) then
                L(i, j) = LU(i, j)
            else
                U(i, j) = LU(i, j)
            end if
        end do
    end do
    ! need to permute B
    do i = 1, ma
        Y(i, :) = B(P(i), :)
    end do
    B = Y
    Y = 0.
    !forward sub for Ly = b
    call forwardsub(L, B, Y, ma, mb)

    !Backward sub for Ux = y
    call backsub(U, Y, X, ma, mb)

  end subroutine LUsolve
    
  subroutine forwardsub(L, B, Y, ma, mb) 

    implicit none

    integer, intent(in) :: ma, mb
    real :: L(:, :), B(:, :), Y(:, :)
    integer :: i, j
    Y = 0. 
     
    Y(1, :) = B(1, :)/L(1, 1)
    do i = 2, ma
        Y(i, :) = (B(i, :) - matmul(L(i, 1:i-1), Y(1:i-1, :)))/L(i, i)
    end do

  end subroutine forwardsub

  subroutine cholesky(A, ma, isSingular, tol)

    implicit none

    integer, intent(in) :: ma
    real, intent(in) :: tol
    logical, intent(out) :: isSingular
    integer :: i, j
    real :: A(:, :)

    !check if A is singular
    !how tho lol?
    
    !clears upper triangular region of A
    do j = 2, ma
        do i = 1, j-1
            A(i, j) = 0.0
        end do  
    end do

 
    A(:, 1) = A(:, 1)/sqrt(A(1, 1))
    do j = 2, ma
        A(j, j) = sqrt(A(j, j) - sum(A(j, 1:j-1)**2))
        A(j+1:ma, j:j) = (A(j+1:ma, j:j) - matmul(A(j+1:ma, 1:j-1), transpose(A(j:j, 1:j-1))))/A(j, j)
    end do
            
  end subroutine cholesky

  subroutine LLsolve(L, B, X, ml, mb)
    !note that B and X must be of the same size

    implicit none

    real, intent(in) :: L(:, :)
    integer, intent(in) :: ml, mb
    real :: B(:, :), X(:, :)

    call forwardsub(L, B, X, ml, mb)
    B = X
    call backsub(transpose(L), B, X, ml, mb)

  end subroutine LLsolve

  subroutine householderQR(A, Rvec, ma, na, isSingular, tol)
    !this routine takes A (some m by n matrix) and returns a matrix with an upper triangular R (diagonal entries of R are
    !stores in rvec and the vectors that form Q on/below the diagonal 

    implicit none

    integer, intent(in) :: ma, na
    real, intent(in) :: tol
    logical :: isSingular
    real :: A(:, :), Rvec(:)
    real, allocatable :: x(:, :)
    real :: norm = 0.0
    integer :: i, j

    Rvec = 0.0

    do i = 1, na    
        allocate(x(ma-i+1, 1))
        x = 0
        !take column below the diagonal 
        x(:, 1) = A(i:ma, i)
        !make vector x = sign(x1)twonorm(x)ihat + x
        call twonorm(x(:, 1), norm)

!        print *, "Norm: ", norm
!        call printmat(x, ma-i+1, 1)

        if(norm > 10.d-15) then
!            print *, "entered if statement"
            x(1, 1) = x(1, 1) + sign(norm, x(1, 1))
            !normalize x
            call twonorm(x(:, 1), norm)
            x = x/norm
        end if
        ! Multiply A by the householder reflector
        A(i:ma, i:na) = A(i:ma, i:na) - 2*matmul(x, matmul(transpose(x), A(i:ma, i:na)))
        Rvec(i) = A(i, i)
        A(i, i) = 0
        A(i:ma, i) = x(:, 1)
        deallocate(x)
        
    end do

  end subroutine householderQR

  subroutine ident(I, ma)
    ! returns an ma x ma identity matrix. 

    implicit none

    integer, intent(in) :: ma
    real :: I(:, :)
    integer :: k

    I = 0.0
    do k = 1, ma
        I(k, k) = 1.0
    end do  

  end subroutine ident

  subroutine formQstar(A, Q, ma, na)
    ! takes A such that on/below diagonal the vectors vi that compose qi are in the ith column of A. 
    ! returns Q = Q1...QN

    implicit none

    integer, intent(in) :: ma, na
    integer :: i, j
    real :: A(:, :),  Q(:, :)
    real, dimension(ma, ma) :: eye, Qi, QTEMP
    
    ! Q = Q1...QN

    ! Qi = | I 0 |
    !      | 0 H |
    
    ! H = I - 2vv^T

    ! QTEMP = 2
    ! from i = n to 1
    ! QTEMP = Qi QTEMP
    
    eye = 0.0
    call ident(eye, ma)
    
    QTEMP = eye

    do j = 1, na
        i = na - j + 1
        Qi = eye
        Qi(i:ma, i:ma) = Qi(i:ma, i:ma) - 2*matmul(A(i:ma, i:i), transpose(A(i:ma, i:i))) ! H = I - 2vv^T
        QTEMP = matmul(Qi, QTEMP)
    end do
    Q = QTEMP(:, 1:na) 

  end subroutine formQstar
    
  subroutine formR(A, R, rvec, na)
    ! takes in A such that above the diagonal are the above diagonal entries of R, and the diagonal entries are in rvec. 
    ! returns and na x na matrix R

    implicit none

    integer, intent(in) :: na
    integer :: i, j
    real :: A(:, :), R(:, :), rvec(:)

    R = 0.0 

    do i = 1, na
        if(i .eq. 1) then
            R(1, 1) = rvec(1)
        else
            !take A above the diagonal
            R(1:i-1, i) = A(1:i-1, i)
            R(i, i) = rvec(i)
        end if
    end do

  end subroutine formR

  subroutine checksingular(A, ma, isSingular, tol)
    ! this hasn't been implemented yet :)
    ! See ya later Nate!
    real :: A(:, :)
    integer, intent(in) :: ma
    logical :: isSingular
    real, intent(in) :: tol

    !checks the diagonal elements to see if they are zero (this would influence the backsub/forward sub
    !a better way would be to 
    
  end subroutine checksingular

  subroutine frobnorm(A, norm)
    !computes ||A|| = sqrt(trace(A^TA))
    implicit none
    real, intent(in) :: A(:, :)
    real, intent(out) :: norm
    norm = sqrt(sum(A**2))
  end subroutine frobnorm

  subroutine tridiagonal(A, ma)

    implicit none

    integer, intent(in) :: ma
    real :: A(:, :)
    real, allocatable :: x(:, :)
    real :: norm
    integer :: i, j


    do i = 1, ma-2
        allocate(x(ma-i, 1))
        x = 0
        x(:, 1) = A(i+1:ma, i)
        call twonorm(x(:, 1), norm)
        x(1, 1) = x(1, 1) + (x(1, 1)/abs(x(1, 1)))*norm
        call twonorm(x(:, 1), norm)
        x = x/norm
        A(i+1:ma, i:ma) = A(i+1:ma, i:ma) - 2*matmul(x, matmul(transpose(x), A(i+1:ma, i:ma)))
        A(1:ma, i+1:ma) = A(1:ma, i+1:ma) - 2*(matmul(A(1:ma, i+1:ma), matmul(x, transpose(x))))
        deallocate(x)
        
    end do

  end subroutine tridiagonal

  subroutine diag(A, ma, D)

    implicit none

    integer :: i
    integer, intent(in) :: ma
    real :: A(:, :), D(:, :)

    D = 0.0

    do i = 1, ma
        D(i, i) = A(i, i)
    end do

  end subroutine diag
  
  subroutine eigQR(A, ma, shift, tol)

    implicit none

    integer, intent(in) :: ma
    real, intent(in) :: tol
    logical, intent(in) :: shift
    real :: A(:, :)
    real, dimension(ma) :: rvec
    real, dimension(ma, ma) :: D, D2
    real, dimension(ma, ma) :: Q, R, eye
    real :: mu, error, norm
    logical :: isSingular = .false.
    integer :: i, j, k
        
    eye = 0.0
    call ident(eye, ma)
    norm = 10000

    if (shift) then

        call tridiagonal(A, ma)    
        i = 0
        do while (norm > tol)
        
            mu = A(ma, ma)

            R = 0.0
            Q = 0.0
            A = A-mu*eye

            call householderQR(A, rvec, ma, ma, isSingular, tol)
            call formR(A, R, rvec, ma)
            call formQstar(A, Q, ma, ma)

            A = matmul(R, Q) + mu*eye

            call diag(A, ma, D2)
            call frobnorm(D2-D, norm)
            D = D2
            i = i + 1
        end do
        print *, "Algorithm converged within "//trim(str(i))//" iterations"

    else

        i = 0
        do while (norm > tol)

            R = 0.0
            Q = 0.0
            call householderQR(A, rvec, ma, ma, isSingular, tol)
            call formR(A, R, rvec, ma)
            call formQstar(A, Q, ma, ma)

            !This routine fails after the matrix becomes diagonal 

            A = matmul(R, Q)

            !compare eigenvals of Ak+1 to Ak
            !note that eye is equal to Dk
            !We use Q's memory for D_k+1
            call diag(A, ma, D)
            call frobnorm(D-D2, norm)
            D2 = D
            i = i + 1
        end do
        print *, "Algorithm converged within "//trim(str(i))//" iterations"
    end if

  end subroutine eigQR

 
  subroutine inviter(A, eig, v, ma, tol)
        
    implicit none

    integer, intent(in) :: ma
    real :: A(:, :), eig, v(:, :), error, tol, mu
    real, dimension(ma, ma) :: As, eye, Q, R
    real, dimension(ma, 1) :: w
    real, dimension(ma) :: rvec
    integer, dimension(ma) :: p
    logical :: dummy
    integer :: i

    v = 0.0


    !This initializes v0 to be a vector of norm 1
    do i = 1, ma
        v(i, 1) = 1.0/sqrt(real(ma))
    end do
        
    call ident(eye, ma)
    error = 10.0
    As = A

    do while (error > tol)

        w = 0.0
   
        mu = eig - tol
        A = As - mu*eye
        
        !Solving the linear system
        call householderQR(A, rvec, ma, ma, dummy, tol)
        call formR(A, R, rvec, ma)
        call formQstar(A, Q, ma, ma)
        v = matmul(transpose(Q), v)
        call backsub(R, v, w, ma, ma)
   
        !using the error memory as a norm funciton for now 
        call twonorm(w(:, 1), error)
        v = w/error

        w = matmul(As, v) - eig*v

        call twonorm(w(:, 1), error)

    end do
 
  end subroutine inviter

end module LinAl
