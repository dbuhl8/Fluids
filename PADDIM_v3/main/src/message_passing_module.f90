MODULE message_passing_module
  USE defprecision_module
  IMPLICIT NONE
  SAVE
!
!*********************************************************************
!***                        GLOBAL DATA                            ***
!*********************************************************************
  INTEGER(kind=kiMPI) :: myid
  INTEGER(kind=kiMPI) :: numtasks
  INTEGER(kind=kiMPI) :: nprocs1,nprocs2    !#tasks of 1st & 2nd transpose

CONTAINS

!*********************************************************************
!***                       Subroutines                             ***
!*********************************************************************

#include "message_passing/start_mpi.f90"
#include "message_passing/stop_mpi.f90"

end MODULE message_passing_module
