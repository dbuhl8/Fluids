#include "defs_MPI.h"
MODULE integral_module
! +--------------------------------------------------------------------------+
! | This module containes routines that apply several integral operators to  |
! | a given scalar field in physical space.                                  |
! | Only volume integration is currently implemented.                        |
! +--------------------------------------------------------------------------+
  USE defprecision_module
  USE parameter_module, ONLY :
  USE mpi_transf_module, ONLY :
  USE message_passing_module, ONLY :
  IMPLICIT NONE
  SAVE
!**************************************************************
!****  Global Data for needed for computing differentials  ****
!**************************************************************
  ! no global data

CONTAINS

!*********************************************
!****  subroutines for taking derivatives ****
!*********************************************

#include "integral/volume_integral.f90"

END MODULE integral_module
