#include "defs_MPI.h"
MODULE testing_module
  USE defprecision_module
  USE state_module, ONLY :
  USE mpi_transf_module, ONLY:
  USE parameter_module, ONLY:
  USE diagnostics_module, ONLY:
  IMPLICIT NONE

CONTAINS

#include "testing/error_Temp_adv_diff.f90"
#include "testing/peak_div_u.f90" 
#include "testing/peak_div_B.f90" 

END MODULE testing_module
