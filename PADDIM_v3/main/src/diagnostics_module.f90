#include "defs_MPI.h"
MODULE diagnostics_module
!+-------------------------------------------------------------------------+
!| This module containes subroutines used for calculating diagnostical     |
!| output parameters (various energies, rms-values, ...).                  |
!|                                                                         |
!| Author: Stephan Stellmach                                               |
!+-------------------------------------------------------------------------+
!
  USE defprecision_module
  USE defs_2D_3D_module  
  USE parameter_module, ONLY :
  USE state_module, ONLY :
  USE integral_module, ONLY :
  USE mpi_transf_module, ONLY :
  IMPLICIT NONE
  SAVE
!****************************************************
!****                   Global Data              ****
!****************************************************
!
! no global data so far
!
!*********************************************************************
!***            Interfaces for generic procedures                  ***
!*********************************************************************
INTERFACE rms
   ! This generic function computes the RMS values of a field, no matter if scalar or vector
   ! or if given in physical or spectral space. 
   MODULE PROCEDURE rms_scalar,rms_vector,rms_spec_scalar,rms_spec_vector
END INTERFACE

CONTAINS

!****************************************************
!****                Subroutines                 ****
!****************************************************
#include "diagnostics/rms_scalar.f90"
#include "diagnostics/rms_spec_scalar.f90"
#include "diagnostics/rms_vector.f90"
#include "diagnostics/rms_spec_vector.f90"
#include "diagnostics/compute_uTCB_rms.f90"
#include "diagnostics/compute_uTCB_minmax.f90"
#include "diagnostics/compute_average_flux.f90"
#include "diagnostics/dissipation_buo.f90"
#include "diagnostics/compute_z_autocorel.f90"

end MODULE diagnostics_module
