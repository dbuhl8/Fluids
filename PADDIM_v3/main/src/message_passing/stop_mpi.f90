SUBROUTINE stop_mpi
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif

  INTEGER(kind=ki) :: ierr
! finalize mpi
  CALL MPI_FINALIZE(ierr)

END SUBROUTINE stop_mpi

