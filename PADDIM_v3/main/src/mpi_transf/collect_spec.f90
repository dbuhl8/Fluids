! collect spectra at process sourcetask for pencil decomposition
! not optimized for efficiency!
SUBROUTINE collect_spec(local,global,sourcetask)
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  COMPLEX(kind=kr)                 :: local(:,:,:)
  COMPLEX(kind=kr)                 :: global(0:,0:,0:)
  INTEGER(kind=kiMPI)              :: sourcetask
  INTEGER(kind=kiMPI),ALLOCATABLE  :: recvcounts(:),rdispl(:)
  INTEGER(kind=ki),ALLOCATABLE     :: sx(:),ex(:),sy(:),ey(:)
  COMPLEX(kind=kr),ALLOCATABLE     :: work(:)
  INTEGER(kind=kiMPI)              :: sendcount,ierr,comm
  INTEGER(kind=ki)                 :: numtasks,myid,process,count
  INTEGER(kind=ki)                 :: i,j,k

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

  ! collect information from other processes in work array
  IF (myid.EQ.sourcetask) THEN     
     ALLOCATE(work((Lmax+1)*MAX((2*Mmax),1)*(2*Nmax)))
     ALLOCATE(recvcounts(0:numtasks-1))
     ALLOCATE(rdispl(0:numtasks-1))
     recvcounts(:) = (2*Nmax) * (ex(:) - sx(:) + 1) * (ey(:) - sy(:) + 1) 
     count = 1 
     DO process = 0,numtasks-1
        rdispl(process) = count - 1
        count = count +  ( ey(process)-sy(process)+1 ) &
        &              * ( ex(process)-sx(process)+1 ) &
        &              * ( 2*Nmax )
     ENDDO
  ELSE
     ALLOCATE(work(1),recvcounts(1),rdispl(1)) ! Some compilers complain in these are not allocated
  ENDIF

  sendcount = (2*Nmax) * (myex_spec - mysx_spec + 1) * (myey_spec - mysy_spec + 1)
  CALL MPI_GATHERV(local,sendcount,PM_MPI_COMPLEX_TYPE,              &
  &                work,recvcounts,rdispl,PM_MPI_COMPLEX_TYPE,       &
  &                sourcetask,comm,ierr )

  ! sort stuff into array global
  IF (myid.EQ.sourcetask) THEN
     count = 1 
     DO process = 0,numtasks-1
        DO j = sy(process),ey(process)
           DO i = sx(process),ex(process)
              DO k = 0,2*Nmax-1
                 global(k,i,j) = work(count)
                 count = count + 1 
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  DEALLOCATE(sx,ex,sy,ey)
  DEALLOCATE(work,recvcounts,rdispl)

END SUBROUTINE collect_spec
