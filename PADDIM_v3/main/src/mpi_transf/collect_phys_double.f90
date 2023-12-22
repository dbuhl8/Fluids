! collect physical fields at process sourcetask for pencil decomposition
! not optimized for efficiency!
SUBROUTINE collect_phys_double(local,global,sourcetask)
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  REAL(kind=krd)                   :: local(:,:,:)
  REAL(kind=krd)                   :: global(0:,0:,0:)
  INTEGER(kind=kiMPI)              :: sourcetask
  INTEGER(kind=kiMPI),ALLOCATABLE  :: recvcounts(:),rdispl(:)
  INTEGER(kind=ki),ALLOCATABLE     :: sy(:),ey(:),sz(:),ez(:)
  REAL(kind=krd),ALLOCATABLE       :: work(:)
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
     ALLOCATE(sy(0:numtasks-1))
     ALLOCATE(ey(0:numtasks-1))
     ALLOCATE(sz(0:numtasks-1))
     ALLOCATE(ez(0:numtasks-1))
  ELSE 
     ALLOCATE(sy(1),ey(1),sz(1),ez(1)) !Some compilers complain if these are not allocated at all
  ENDIF
     

  CALL MPI_GATHER(mysy_phys,1,MPI_INTEGER,sy,1,MPI_INTEGER,              &
  &               sourcetask,comm,ierr)
  CALL MPI_GATHER(myey_phys,1,MPI_INTEGER,ey,1,MPI_INTEGER,              &
  &               sourcetask,comm,ierr)
  CALL MPI_GATHER(mysz_phys,1,MPI_INTEGER,sz,1,MPI_INTEGER,              &
  &               sourcetask,comm,ierr)
  CALL MPI_GATHER(myez_phys,1,MPI_INTEGER,ez,1,MPI_INTEGER,              &
  &               sourcetask,comm,ierr)

  ! collect information from other processes in work array
  IF (myid.EQ.sourcetask) THEN     
     ALLOCATE(work(Nx*Ny*Nz))
     ALLOCATE(recvcounts(0:numtasks-1))
     ALLOCATE(rdispl(0:numtasks-1))
     recvcounts(:) = Nx * (ey(:) - sy(:) + 1) * (ez(:) - sz(:) + 1) 
     count = 1 
     DO process = 0,numtasks-1
        rdispl(process) = count - 1
        count = count +  ( ez(process)-sz(process)+1 ) &
        &              * ( ey(process)-sy(process)+1 ) &
        &              * Nx
     ENDDO
  ELSE
     ALLOCATE(work(1),recvcounts(1),rdispl(1)) !Some compiler coamplain if these are not allocated at all
  ENDIF


  sendcount = Nx * (myey_phys - mysy_phys + 1) * (myez_phys - mysz_phys + 1)
  CALL MPI_GATHERV(local,sendcount,PM_MPI_DOUBLE_FLOAT_TYPE,               &
  &                work,recvcounts,rdispl,PM_MPI_DOUBLE_FLOAT_TYPE,        &
  &                sourcetask,comm,ierr )

  ! sort stuff into array global
  IF (myid.EQ.sourcetask) THEN
     count = 1 
     DO process = 0,numtasks-1
        DO k = sz(process),ez(process)
           DO j = sy(process),ey(process)
              DO i = 0,Nx-1
                 global(i,j,k) = work(count)
                 count = count + 1 
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  DEALLOCATE(sy,ey,sz,ez)
  DEALLOCATE(work,recvcounts,rdispl)

END SUBROUTINE collect_phys_double

