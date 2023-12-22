SUBROUTINE init_transpose_1d_decomp(info,comm,N1,N2,M3)
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  INTEGER(kind=ki)               :: N1,N2,M3
  INTEGER(kind=kiMPI)            :: comm
  TYPE(transpose_info_1d_decomp) :: info
  INTEGER(kind=ki)               :: ierr

  ! dimensions of the array which is to be transposed 
  info%N1 = N1
  info%N2 = N2
  info%M3 = M3
  ! set communicator assiciated with this transpose
  info%communicator = comm
  ! find out my id in communicator comm
  CALL MPI_COMM_RANK(comm, info%myid, ierr)
  ! how many others are there in the group associated with communicator comm?
  CALL MPI_COMM_SIZE(comm, info%numtasks, ierr)
  ! allocate memory to store subdomain information
  ALLOCATE(info%s1(0:info%numtasks-1))
  ALLOCATE(info%e1(0:info%numtasks-1))
  ALLOCATE(info%s2(0:info%numtasks-1))
  ALLOCATE(info%e2(0:info%numtasks-1))
  ! compute decomposion of first dimension among numtask tasks
  CALL decomp_1d(N1,info%numtasks,info%myid,info%mys1,info%mye1)
  ! compute decomposion of second dimension among numtask tasks
  CALL decomp_1d(N2,info%numtasks,info%myid,info%mys2,info%mye2)
  ! find out how large the arrays of the other processors are and store in info
  CALL MPI_ALLGATHER(info%mys1,1,MPI_INTEGER,info%s1,1,MPI_INTEGER,comm,ierr)
  CALL MPI_ALLGATHER(info%mye1,1,MPI_INTEGER,info%e1,1,MPI_INTEGER,comm,ierr)
  CALL MPI_ALLGATHER(info%mys2,1,MPI_INTEGER,info%s2,1,MPI_INTEGER,comm,ierr)
  CALL MPI_ALLGATHER(info%mye2,1,MPI_INTEGER,info%e2,1,MPI_INTEGER,comm,ierr)
  ! allocate memory to store integer arrays passed to MPI_ALL_TO_ALL
  ALLOCATE(info%sendcounts(0:info%numtasks-1))
  ALLOCATE(info%sdispl(0:info%numtasks-1))
  ALLOCATE(info%recvcounts(0:info%numtasks-1))
  ALLOCATE(info%rdispl(0:info%numtasks-1))
  ! initialize these integer arrays
  info%sendcounts(:) =   info%M3 * (info%e1(:)-info%s1(:)+1)               &
  &                              * (info%mye2-info%mys2+1)
  info%sdispl(:) = info%M3 * (info%s1(:)-1)                                &
  &                       *  (info%mye2-info%mys2+1)
  info%recvcounts(:) = info%M3 * (info%e2(:)-info%s2(:)+1)                 &
  &                            * (info%mye1-info%mys1+1)
  info%rdispl(:) = info%M3 * (info%s2(:) - 1)                              &
  &                        * (info%mye1-info%mys1+1)

END SUBROUTINE init_transpose_1d_decomp
