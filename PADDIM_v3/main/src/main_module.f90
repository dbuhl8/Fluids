#include "defs_MPI.h"
module main_module
!+--------------------------------------------------------------------+
!| This module contains the subroutine needed by the main part of     |
!| the program.                                                       |
!+--------------------------------------------------------------------+
  use defprecision_module
  USE defs_2D_3D_module
  USE parameter_module, ONLY :
  USE mpi_transf_module, ONLY :
  USE pnetCDF_IO_module, ONLY :
  USE IO_module, ONLY :
  USE state_module, ONLY :
  use message_passing_module, ONLY :
  use diff_op_module, ONLY :
  IMPLICIT NONE
  SAVE
!*********************************************************************
!***                        GLOBAL DATA                            ***
!*********************************************************************

  ! Coefficients in AB/BDF3 time stepping scheme
  real(kind=kr) :: ta0,ta1,ta2,ta3,tb1,tb2,tb3 

CONTAINS

!*********************************************************************
!***                       Subroutines                             ***
!*********************************************************************

! 2nd order Runge-Kutta / Crank-Nicholson (Startup scheme for AB/BDF3)
#include "main/RK2_CN/timestep_RK2_CN.f90"
#include "main/RK2_CN/step1_RK2_CN.f90"
#include "main/RK2_CN/step2_RK2_CN.f90"
! third order Adams-Bashforth / Backward-Differencing (semi-implicit multi-step)
#include "main/AB_BDF3/timestep_AB_BDF3.f90"
#include "main/AB_BDF3/tmstp_buoyancy_AB_BDF3.f90"
#include "main/AB_BDF3/tmstp_velocity_AB_BDF3.f90"
#include "main/AB_BDF3/tmstp_magnetic_AB_BDF3.f90"
#include "main/AB_BDF3/comp_coeff_AB_BDF3.f90"
#include "main/AB_BDF3/CFL_AB_BDF3.f90"
#include "main/AB_BDF3/adapt_dt_AB_BDF3.f90"
! universal routines 
#include "main/read_parameter.f90"
#include "main/init.f90"
#include "main/free_allocated_memory.f90"
#include "main/crhs_buoyancy.f90"
#include "main/crhs_velocity.f90"
#include "main/crhs_magnetic.f90"
#include "main/compute_phys_space_vars.f90"


end module main_module
