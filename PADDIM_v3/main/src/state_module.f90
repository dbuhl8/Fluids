!+--------------------------------------------------------------------+
!| This module contains datatypes and global data needed to store     |
!| the variable which determine the physical state of the system.     |
!| Data from previous timesteps (old nonlinear terms,...) which are   |
!| needed for the time stepping are also included.                    |
!| Furthermore, some subroutines are included which perform data      |
!| initializations, etc.                                              |
!+--------------------------------------------------------------------+
MODULE state_module
  USE defprecision_module
  USE defs_2D_3D_module
  USE parameter_module, ONLY :
  USE mpi_transf_module, ONLY :
  USE diff_op_module, ONLY : 
  IMPLICIT NONE
  SAVE
  PUBLIC :: velocity,buoyancy
!*********************************************************************
!***                      DERIVED DATATYPES                        ***
!*********************************************************************
  TYPE velocity 
    REAL(kind=kr) ,   DIMENSION(:,:,:,:), POINTER :: phys
    COMPLEX(kind=kr), DIMENSION(:,:,:,:,:), POINTER :: spec
    REAL(kind=kr) ,   DIMENSION(:,:,:,:), POINTER :: curl
    COMPLEX(kind=kr), DIMENSION(:,:,:,:,:), POINTER :: rhs
  END TYPE velocity

  ! PH{
  TYPE field 
    REAL(kind=kr) ,   DIMENSION(:,:,:,:), POINTER :: phys
    COMPLEX(kind=kr), DIMENSION(:,:,:,:,:), POINTER :: spec
    REAL(kind=kr) ,   DIMENSION(:,:,:,:), POINTER :: curl
    COMPLEX(kind=kr), DIMENSION(:,:,:,:,:), POINTER :: rhs
  END TYPE field
  ! }PH
  
  TYPE buoyancy
    REAL(kind=kr) ,   DIMENSION(:,:,:), POINTER :: phys
    COMPLEX(kind=kr), DIMENSION(:,:,:,:), POINTER :: spec
    COMPLEX(kind=kr), DIMENSION(:,:,:,:), POINTER :: rhs
  END TYPE buoyancy

!*********************************************************************
!***                        GLOBAL DATA                            ***
!*********************************************************************
  INTEGER(kind=ki),PARAMETER :: time_levels = 3        ! Number of time levels that 
                                                       ! need to be stored
! In the spec arrays, the spectrum is saved at multiple time levels
! The last index specefies the time level:
  INTEGER (kind=ki) :: ltime0 = 1 !spectrum at time step t_n
  INTEGER (kind=ki) :: ltime1 = 2 !spectrum at time step t_{n-1}
  INTEGER (kind=ki) :: ltime2 = 3 !spectrum at time step t_{n-2}

  INTEGER (kind=ki) :: rtime1 = 1 !RHS at time step t_{n-1}
  INTEGER (kind=ki) :: rtime2 = 2 !RHS at time step t_{n-2}

CONTAINS

#include "state/allocate_uTCB.f90"
#include "state/deallocate_uTCB.f90"
#include "state/set_initial_condition.f90"
#include "state/init_u_phys.f90"
#include "state/init_B_phys.f90"
#include "state/init_Temp_phys.f90"
#include "state/init_Chem_phys.f90"
#include "state/shift_time_pointers.f90"
#include "state/decomp_independent_random.f90"
#include "state/make_solenoidal.f90"
END MODULE state_module
