SUBROUTINE create_pencil_MPI_communicators(nprocs1,nprocs2,comm1,comm2,basic_communicator)
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  INTEGER(kind=ki)    :: nprocs1,nprocs2
  INTEGER(kind=kiMPI) :: comm1,comm2,basic_communicator
  INTEGER(kind=ki)    :: myid,numtasks
  INTEGER(kind=ki)    :: ierr
  INTEGER(kind=ki)    :: color1,color2,key1,key2

  CALL MPI_COMM_SIZE(basic_communicator, numtasks, ierr)
  CALL MPI_COMM_RANK(basic_communicator, myid, ierr)

  ! check in nprocs2 * nprocs2 = total number of processors
  IF (nprocs1*nprocs2 .NE. numtasks) THEN
     PRINT*,"Error: Total number of tasks different from nprocs1 * nprocs2 !"
     CALL MPI_FINALIZE(ierr)
     STOP
  ENDIF

  ! Define two communicators. The associated process groups contain
  ! processes witch are in the same xy process plane or in the same 
  ! yz process plane. Using the new communicators, each transpose 
  ! can be reduced to 1d transposes.

  color1 = myid / nprocs1
  key1   = MOD(myid,nprocs1)

  color2 = MOD(myid,nprocs1)
  key2   = myid / nprocs1

  CALL MPI_COMM_SPLIT(basic_communicator,color1,key1,comm1,ierr)
  CALL MPI_COMM_SPLIT(basic_communicator,color2,key2,comm2,ierr)

END SUBROUTINE create_pencil_MPI_communicators
