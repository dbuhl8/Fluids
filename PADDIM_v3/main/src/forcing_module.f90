!+--------------------------------------------------------------------+
!| This module is responsible for computing the forcing for PADDI     |
!| force_spec is ultimatelty returned to crhs_velocity.f90            |
!| This module regularly uses kx, ky, kz and only shares              |
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


real(kind=kr), allocatable :: gpTimeVals(:), gpForcingVals(:, :), interpSlopeVals(:, :)
real(kind=kr) :: gaussian_tmscl, tol, tstep
integer(kind=ki), allocatable :: waveNumMap(:, :)
integer(kind=ki) :: cTimeIndex, numGPcolumns, numGProws, window_pts, window_skip
integer(kind=ki) :: numTrashLines
real(kind=kr)    :: KMAX_forcing
character*20 :: dataFormat
character*20 :: headerFormatReal, headerFormatInt


CONTAINS

#include "forcing/forcing.f90" ! this needs to be fixed later
#include "forcing/gpinterp.f90"
#include "forcing/read_stochastic_forcing.f90"
#include "forcing/close_stochastic_forcing.f90"
#include "forcing/gaussian.f90"

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

END MODULE forcing_module
