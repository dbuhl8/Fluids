



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

integer :: n, n2, i, j
real(8) :: xstart, xstop, xstep, kscale
real(8), allocatable :: xtrain(:), ftrain(:), xsample(:), fsample(:), fact(:)
character(len = 11) :: tab

n = 15
n2 = 100
kscale = sqrt(1.0)

allocate(xtrain(n), ftrain(n), xsample(n2), fsample(n2), fact(n2))

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
    xsample(i) = xstart + (i-1) * xstep
end do
fact = sin(xsample)

!Calling Gaussian Regression Function from gaussian_mod
fsample = dgr(xtrain, ftrain, xsample, n, n2, kscale, 'S', 'E')

!Generating data column tab length for gnuplot reading 
tab = char(11)

!Writing to .dat files
open(10, file = '../plotData/regressiondata.dat')
    do i = 1, n2
        write(10, *) xsample(i),tab, fsample(i),tab, fact(i)
    end do
close(10)

!Two files are needed since the amount of data in their respective vectors is distinct
open(11, file = 'plotData/traindata.dat')
    do i = 1, n
        write(11, *) xtrain(i),tab, ftrain(i)
    end do
close(11)


end program gaussian



