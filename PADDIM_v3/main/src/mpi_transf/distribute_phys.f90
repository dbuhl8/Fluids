! distribute physical fields from process sourcetask for pencil decomposition
! not optimized for efficiency!
SUBROUTINE distribute_phys(global,local,sourcetask)
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  REAL(kind=kr)                   :: local(:,:,:)
  REAL(kind=kr)                   :: global(0:,0:,0:)
  INTEGER(kind=kiMPI)             :: sourcetask
  INTEGER(kind=kiMPI),ALLOCATABLE :: sendcounts(:),sdispl(:)
  INTEGER(kind=ki),ALLOCATABLE    :: sy(:),ey(:),sz(:),ez(:)
  REAL(kind=kr),ALLOCATABLE       :: work(:)
  INTEGER(kind=kiMPI)             :: recvcount,ierr,comm
  INTEGER(kind=ki)                :: numtasks,myid,process,count
  INTEGER(kind=ki)                :: i,j,k

#ifdef TWO_DIMENSIONAL
  comm = two_dim_transp_info%communicator
#else
  comm = pencil_transpose_info_zx_yz%comm_global
#endif

  CALL MPI_COMM_SIZE(comm,numtasks,ierr)
  CALL MPI_COMM_RANK(comm,myid,ierr)

  ! find out what data the other processes hold
  IF (myid.EQ.sourcetask) THEN
     ALLOCATE(sy(0:numtasks-1))
     ALLOCATE(ey(0:numtasks-1))
     ALLOCATE(sz(0:numtasks-1))
     ALLOCATE(ez(0:numtasks-1))
  ELSE
     ALLOCATE(sy(1),ey(1),sz(1),ez(1)) !Some compiler coamplain if these are not allocated at all
  endif

  CALL MPI_GATHER(mysy_phys,1,MPI_INTEGER,sy,1,MPI_INTEGER,              &
  &               sourcetask,comm,ierr      )
  CALL MPI_GATHER(myey_phys,1,MPI_INTEGER,ey,1,MPI_INTEGER,              &
  &               sourcetask,comm,ierr      )
  CALL MPI_GATHER(mysz_phys,1,MPI_INTEGER,sz,1,MPI_INTEGER,              &
  &               sourcetask,comm,ierr      )
  CALL MPI_GATHER(myez_phys,1,MPI_INTEGER,ez,1,MPI_INTEGER,              &
  &               sourcetask,comm,ierr      )

  ! write send packets into continous array work
  IF (myid.EQ.sourcetask) THEN
     ALLOCATE(work(Nx*Ny*(Nz+1)))
     ALLOCATE(sendcounts(0:numtasks-1))
     ALLOCATE(sdispl(0:numtasks-1))

     count = 1 
     DO process = 0,numtasks-1
        sdispl(process) = count - 1
        DO k = sz(process),ez(process)
           DO j = sy(process),ey(process)
              DO i = 0,Nx-1
                 work(count) = global(i,j,k)
                 count = count + 1 
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     sendcounts(:) = Nx * (ey(:) - sy(:) + 1) * (ez(:) - sz(:) + 1) 
  ELSE
     ALLOCATE(work(1),sendcounts(1),sdispl(1)) ! Some compilers complain in these are not allocated
  endif

  ! scatter the data among processes 
  recvcount = Nx * (myey_phys - mysy_phys + 1) * (myez_phys - mysz_phys + 1)
  CALL MPI_SCATTERV(work,sendcounts,sdispl,PM_MPI_FLOAT_TYPE,  &
  &                 local ,recvcount,PM_MPI_FLOAT_TYPE,        &
  &                 sourcetask,comm, &
  &                 ierr )
  
  ! deallocate temporary memory
  DEALLOCATE(sy,ey,sz,ez)
  DEALLOCATE(work)
  DEALLOCATE(sendcounts,sdispl)

END SUBROUTINE distribute_phys
