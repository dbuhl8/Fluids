! distribute spectra from process sourcetask for pencil decomposition
! not optimized for efficiency!
SUBROUTINE distribute_spec(global,local,sourcetask)
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  COMPLEX(kind=kr)                :: local(:,:,:)
  COMPLEX(kind=kr)                :: global(0:,0:,0:)
  INTEGER(kind=kiMPI)             :: sourcetask
  INTEGER(kind=kiMPI),ALLOCATABLE :: sendcounts(:),sdispl(:)
  INTEGER(kind=ki),ALLOCATABLE    :: sx(:),ex(:),sy(:),ey(:)
  COMPLEX(kind=kr),ALLOCATABLE    :: work(:)
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
     ALLOCATE(sx(0:numtasks-1))
     ALLOCATE(ex(0:numtasks-1))
     ALLOCATE(sy(0:numtasks-1))
     ALLOCATE(ey(0:numtasks-1))
  ELSE
     ALLOCATE(sx(1),ex(1),sy(1),ey(1)) !Some compiler coamplain if these are not allocated at all
  ENDIF

  CALL MPI_GATHER(mysx_spec,1,MPI_INTEGER,sx,1,MPI_INTEGER,              &
  &               sourcetask,comm,ierr      )
  CALL MPI_GATHER(myex_spec,1,MPI_INTEGER,ex,1,MPI_INTEGER,              &
  &               sourcetask,comm,ierr      )
  CALL MPI_GATHER(mysy_spec,1,MPI_INTEGER,sy,1,MPI_INTEGER,              &
  &               sourcetask,comm,ierr      )
  CALL MPI_GATHER(myey_spec,1,MPI_INTEGER,ey,1,MPI_INTEGER,              &
  &               sourcetask,comm,ierr      )

  ! write send packets into continous array work
  IF (myid.EQ.sourcetask) THEN
     ALLOCATE(work( (2*Nmax)*(Lmax+1)*MAX((2*Mmax),1) ))
     ALLOCATE(sendcounts(0:numtasks-1))
     ALLOCATE(sdispl(0:numtasks-1))

     count = 1 
     DO process = 0,numtasks-1
        sdispl(process) = count - 1
        DO j = sy(process),ey(process)
           PRINT*,myid,process,j
           DO i = sx(process),ex(process)
              DO k = 0,2*Nmax-1
                 work(count) = global(k,i,j)
                 count = count + 1 
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     sendcounts(:) = (2*Nmax) * (ey(:) - sy(:) + 1) * (ex(:) - sx(:) + 1) 
  ELSE
     ALLOCATE(work(1),sendcounts(1),sdispl(1)) ! Some compilers complain in these are not allocated
  ENDIF

  ! scatter the data among processes 
  recvcount = (2*Nmax) * (myey_spec - mysy_spec + 1) * (myex_spec - mysx_spec + 1)  
  CALL MPI_SCATTERV(work,sendcounts,sdispl,PM_MPI_COMPLEX_TYPE,  &
  &                 local ,recvcount,PM_MPI_COMPLEX_TYPE,        &
  &                 sourcetask,comm, &
  &                 ierr )

  ! deallocate temporary memory
  DEALLOCATE(sx,ex,sy,ey)
  DEALLOCATE(work)
  DEALLOCATE(sendcounts,sdispl)


END SUBROUTINE distribute_spec
