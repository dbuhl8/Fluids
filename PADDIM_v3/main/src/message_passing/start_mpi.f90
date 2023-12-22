SUBROUTINE start_mpi
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  INTEGER(kind=ki) :: ierr

! initialize mpi
  CALL MPI_INIT(ierr)
! find out my id
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
! how many others are there?
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
! wait until all processes now who they are and how many others are there
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE start_mpi
