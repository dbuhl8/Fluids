!+--------------------------------------------------------------------+
!| This module contains functions to apply several differential       |
!| operators to field given in spectral space.                        |
!+--------------------------------------------------------------------+
MODULE diff_op_module
  USE defprecision_module
  USE defs_2D_3D_module
  USE mpi_transf_module, ONLY:
  USE parameter_module, ONLY:
  IMPLICIT NONE
  SAVE
  !*********************************************************************
  !***                        GLOBAL DATA                            ***
  !*********************************************************************
  
  ! no global data
  
CONTAINS

  !*********************************************************************
  !***                       Subroutines                             ***
  !*********************************************************************
  
#include "diff_op/div.f90"
#include "diff_op/curl.f90"
#include "diff_op/d_by_dx.f90"
#include "diff_op/d_by_dy.f90"
#include "diff_op/d_by_dz.f90"

END MODULE diff_op_module
