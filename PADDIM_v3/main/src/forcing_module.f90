!+--------------------------------------------------------------------+
!| This module is responsible for computing the forcing for PADDI     |
!| force_spec is ultimatelty returned to crhs_velocity.f90            |
!+--------------------------------------------------------------------+
#include "defs_MPI.h"
MODULE forcing_module
  USE defprecision_module
  USE defs_2D_3D_module
  USE parameter_module
  USE mpi_transf_module
  USE diff_op_module
  USE MPI
  USE gaussian_mod
  USE message_passing_module, only : myid
    
  IMPLICIT NONE
  SAVE

complex(kind=kr), pointer :: force_spec(:,:,:,:)
real(kind=kr), pointer :: force_real(:,:,:,:)
real(kind=kr), allocatable :: gpTimeVals(:), 
real(kind=kr), allocatable :: interpSlopeVals(:, :)
real(kind=kr), allocatable :: gpForcingVals(:, :) 
real(kind=kr) :: gaussian_tmscl
real(kind=kr) :: tstep
real(kind=kr) :: tol
real(kind=kr) :: KMAX_forcing
integer(kind=ki), allocatable :: waveNumMap(:, :)
integer(kind=ki) :: cTimeIndex
integer(kind=ki) :: numGPcolumns
integer(kind=ki) :: numGProws
integer(kind=ki) :: window_pts
integer(kind=ki) :: window_skip
integer(kind=ki) :: tot_num_forced
integer(kind=ki) :: numTrashLines
character*20 :: dataFormat
character*20 :: headerFormatReal, headerFormatInt
logical :: c2r, forcingCPU, usepf



CONTAINS

#include "forcing/forcing.f90" 
#include "forcing/gpinterp.f90"
#include "forcing/read_stochastic_forcing.f90"
#include "forcing/close_forcing.f90"
#include "forcing/init_forcing.f90"
#include "forcing/gaussian.f90"

!this is a small function I got from stackoverflow to turn integers into strings in fortran. CAREFUL FOR naming conflicts - DB
character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

END MODULE forcing_module
