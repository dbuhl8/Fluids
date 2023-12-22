MODULE parameter_module
!+------------------------------------------------------------------------+
!| This module contains all the parameters necessary for the calculation. |
!| The code solves the following equations:                               |
!|                                                                        |
!| du                                                                     |
!| -- = - \nabla P + D_visc \nabla^2 u + B_therm T e_z - B_comp C e_z     |
!| dt                                                                     |
!|                                                                        |
!| dT                                                                     |
!| -- = D_therm \nabla^2 T + S_therm w                                    |
!| dt                                                                     |
!|                                                                        |
!| dC                                                                     |
!| -- = D_comp \nabla^2 C + S_comp w                                      |
!| dt                                                                     |
!|                                                                        |
!| \nabla * u = 0                                                         |
!|                                                                        |
!| in a triply periodic cube.                                             |
!|                                                                        |
!| Coefficients:                                                          |
!|                                                                        |
!|   B_therm, B_comp represent the coefficients in fromt of the buoyancy  |
!|                   terms;                                               |
!|   C_Lorentz represents the coefficient in front of the Lorentz force   |
!|                   term in  the momentum equation                       |
!|   D_visc,D_therm,D_comp, D_mag are the viscous, thermal, compositional,|
!|                   and magnetic diffusion parameters;                   |
!|   S_therm, S_comp are the coefficients resulting from the background   |
!|                   stratification.                                      | 
!|   R, Theta are the coefficients resulting from the background rotation |
!|                                                                        |
!|                                                                        |
!| Implementing the algorithm for these general equations has the         |
!| advantage that switching between different non-dimensionalizations     |
!| becomes a trivial matter.                                              |  
!+------------------------------------------------------------------------+
  USE defprecision_module
  IMPLICIT NONE 
  SAVE
!*********************************************************************
!***                        GLOBAL DATA                            ***
!*********************************************************************
!
! Parameters which result from user settings
! ------------------------------------------ 
  INTEGER(kind=ki) :: Lmax       ! highest order in x-direction
  INTEGER(kind=ki) :: Mmax       ! highest order in y-direction
  INTEGER(kind=ki) :: Nmax       ! highest irder in Chebyshev expansion
  INTEGER(kind=ki) :: Nx         ! = Number of grid points in x-direction
  INTEGER(kind=ki) :: Ny         ! = Number of grid points in y-direction
  INTEGER(kind=ki) :: Nz         ! = Number of grid points in z-direction
  INTEGER(kind=ki) :: Nsteps     ! Number of time steps 
  REAL(kind=kr)    :: B_therm    ! Thermal buoyancy coefficient
  REAL(kind=kr)    :: B_comp     ! Compositional buoyancy coefficient
  REAL(kind=kr)    :: C_Lorentz  ! Lorentz force coefficient (!PH)
  REAL(kind=kr)    :: D_visc     ! Viscous diffusion parameter
  REAL(kind=kr)    :: D_therm    ! Thermal diffusion parameter
  REAL(kind=kr)    :: D_comp     ! Compositional diffusion parameter
  REAL(kind=kr)    :: D_mag      ! Magnetic diffusion parameter (!PH)
  REAL(kind=kr)    :: S_therm    ! Thermal background stratification parameter
  REAL(kind=kr)    :: S_comp     ! Compositional background stratification parameter
  REAL(kind=kr)    :: R          ! Rotational Parameter
  REAL(kind=kr)    :: Theta      ! Angle between rotational axis and gravity
  REAL(kind=kr)    :: Gammax     ! length (x-direction) of the box
  REAL(kind=kr)    :: Gammay     ! width (y-direction) of the box
  REAL(kind=kr)    :: Gammaz     ! height (z-direction) of the box
  REAL(kind=kr)    :: dt_initial ! Initial time step
  REAL(kind=kr)    :: dt_max     ! maximum time step
  REAL(kind=kr)    :: CFL_safety_fac ! Safety Factor for CFL-Condition
  REAL(kind=kr)    :: alpha      ! =2 * pi / Gammax
  REAL(kind=kr)    :: beta       ! =2 * pi / Gammay
  REAL(kind=kr)    :: gamma      ! =2 * pi / Gammaz
  ! FFTW wisdom stuff
  CHARACTER(LEN=100):: FFTW_wisdom_file     ! File to store FFTW wisdom
  LOGICAL           :: use_FFTW_wisdom_file ! Should wisdom be imported from this file?
  INTEGER(kind=ki),PARAMETER :: uwisdom=60  ! File unit of FFTW widsom file
  ! FFT storage scheme
  REAL (kind=kr) , POINTER   :: kx(:),ky(:),kz(:) !Mapping tables - array index to wavenumber kx(i) = i * alpha, ky, kz

!
! Constants which are frequently needed by the program
! ----------------------------------------------------

  REAL   (kind=kr),PARAMETER :: pi = 4._kr*ATAN(1._kr) ! Pi

CONTAINS
!*********************************************************************
!***                       Subroutines                             ***
!*********************************************************************
  
#include "parameter/allocate_fft_storage_scheme.f90"
#include "parameter/init_fft_storage_scheme.f90"
#include "parameter/deallocate_fft_storage_scheme.f90"

END MODULE parameter_module

