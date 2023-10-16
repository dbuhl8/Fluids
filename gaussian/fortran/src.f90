


!Fortran File

!Written by Dante Buhl

!The purpose of this file is to create a gaussian random process and display a graph using gnuplot

! This code will be using LAPACK in order to perform an eigen decomposition in order to compute the covariance matrix necessary for
! to generate a multivariate gaussian distribution. The goal is to accurately use the DSYEV function from lapack and to make a
! gnuplot which models a gaussian random process. 

! I need to relearn some fortran in order to make this happen. I do however know what my code should look like.

program gaussian

use gaussian_mod

implicit none

integer :: n, n2

n = 15
n2 = 100

print *, gausthing(n, n2)

contains 

function gausthing(n, n2) result(l)
implicit none
!variables
integer, intent(in) :: n, n2

integer :: i, j, l
real(8) :: xstart, xstop, xstep
real(8), dimension(n) :: xtrain, ftrain
real(8), dimension(n2) :: mu_star, ftest, fact, xtest
real(8), dimension(n, n) :: sigma, inv_sigma
real(8), dimension(n, n2) :: sigmak
real(8), dimension(n2, n) :: sigmakt
real(8), dimension(n2, n2) :: sigmakk, sigmastar, sigma_fact
character(len = 11) :: tab
integer, external :: ilaenv


!Allocating memory for our matrices
!allocate(xtrain(n), ftrain(n), xtest(n2), ftest(n2), fact(n2))
!allocate(mu_star(n2), sigma(n, n), sigmak(n, n2), sigmakt(n2, n), sigmakk(n2, n2))
!allocate(inv_sigma(n, n), sigmastar(n2, n2), sigma_fact(n2, n2))

l = n
!Generating Random Samples for xtrain, xtest
call random_number(xtrain)

!Fitting random sample to the desired domain
xtrain = xtrain * 2 * pi
ftrain = sin(xtrain)

!Generating Values for xtest
xstart = 0
xstop = 2 * pi
xstep = (xstop - xstart)/n2
do i = 1, n2
    xtest(i) = xstart + (i-1) * xstep
end do
fact = sin(xtest)



!Generating Covariance Matrices with our Kernal
sigma = esqk(xtrain, xtrain, n, n)
sigmak = esqk(xtrain, xtest, n, n2)
sigmakt = transpose(sigmak)
sigmakk = esqk(xtest, xtest, n2, n2)


!Generating INV_SIGMA using eigendecompision
!inv_sigma = sigma
!call dsyev('V', 'L', n, inv_sigma, n, inv_eigen, work1, lwork1, info)
!cond_fact = inv_eigen(n) * 10**(1)
!do i = 1, n
!    if(inv_eigen(i) .le. cond_fact) then
!        diag2(i, i) = 0
!    else 
!        diag2(i, i) = 1/inv_eigen(i)
!    end if
!end do 
!inv_sigma = matmul(inv_sigma, matmul(diag2, transpose(inv_sigma)))

inv_sigma = isqla(sigma, n)

!### THIS INVERSE PROCEDURE LOOKS GOOD the issue isn't the inversion ###
!Check to see if inverse is computed correctly
!sigma = matmul(inv_sigma, sigma)
!open(15, file='invCheck.dat')
!do i = 1, n
!    write(15, *) sigma(i, :)
!end do
!close(15)
!The issue is either in the kernal function or in the ftrain issue
!The kernal will be rewritten to see if it is the issue
!Computing the Mu vector using matrix multiplication. Have to compute the inverse first
mu_star = matmul(sigmakt, matmul(inv_sigma, ftrain))

!Computing Covariance Matrix for our test sample
sigmastar = sigmakk - matmul(sigmakt, matmul(inv_sigma, sigmak))

!Eigen decomposition of covariance matrix 
!sigma_fact = sigmastar

!V becomes the matrix of eigenvectors, W is an vector of eigenvalues
!call dsyev('V', 'L', n2, sigma_fact, n2, W, work2, lwork2, info)

!Creating the diagonal eigenvalue matrix
!diag = 0.0
!do i = 1, n2
!    diag(i, i) = sqrt(w(i))
!end do

!Creating the "square-root" matrix of sigmastar
!sigma_fact = matmul(sigma_fact, diag)

!Getting rid of NAN values in the matrix
!do j = 1, n2
!    do i = 1, n2
!        !Variables assigned to NAN can't even equal themselves
!        if (sigma_fact(i, j) .ne. sigma_fact(i, j)) then
!            sigma_fact(i, j) = 0
!        end if
!    end do
!end do

sigma_fact = sysq(sigmastar, n2)

!Generating Gaussian Samples
ftest = rnorm(n2)
ftest = mu_star + matmul(sigma_fact, ftest)
!You can generate more samples with more calls to rnorm

!Generating data column tab length for gnuplot reading 
tab = char(11)

!Writing to .dat files
open(10, file = '../plotData/regressiondata.dat')
    do i = 1, n2
        write(10, *) xtest(i),tab, mu_star(i),tab, fact(i),tab, ftest(i)
    end do
close(10)

!Two files are needed since the amount of data in their respective vectors is distinct
open(11, file = '../plotData/traindata.dat')
    do i = 1, n
        write(11, *) xtrain(i),tab, ftrain(i)
    end do
close(11)

end function gausthing


end program gaussian



