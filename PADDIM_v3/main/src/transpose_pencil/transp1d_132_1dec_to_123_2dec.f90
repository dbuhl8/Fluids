! The follwing subroutine transposes a three dimensional 
! array in(:,:,:) which is decomposed along the x1-direction to an 
! array out(:,:,:) which is decomposed along the x2-direction. 
! The order in which the different directions are stored changes 
! from x1:x3:x2 -> x1:x2:x3.  
SUBROUTINE transp1d_132_1dec_to_123_2dec(in,out,info)
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  TYPE(transpose_info_1d_decomp) :: info
  COMPLEX(kind=kr)               :: in(info%mys1:,:,:)
  COMPLEX(kind=kr)               :: out(:,info%mys2:,:)
  COMPLEX(kind=kr),ALLOCATABLE   :: work(:)
  INTEGER(kind=ki)               :: ierr
  INTEGER(kind=ki)               :: j,k,cell
  INTEGER(kind=kli)              :: count,work_length
! allocate temporary memory 
! -------------------------
  work_length =  int(info%mye2-info%mys2+1,kind=kli)  &
  &             *int(info%N1,kind=kli)                &
  &             *int(info%M3,kind=kli)
  ALLOCATE(work(work_length))
! perform the all to all communication using MPI_ALLTOALLV
! --------------------------------------------------------
  CALL MPI_ALLTOALLV(in  ,info%recvcounts,info%rdispl,PM_MPI_COMPLEX_TYPE,  &
  &                  work,info%sendcounts,info%sdispl,PM_MPI_COMPLEX_TYPE,  &
  &                  info%communicator,ierr)
  ! note that recvcounts und rdispl for the backward transposition 
  ! are the same as sendcounts and senddispl for the forward transposition
  ! and vice versa
! sort received packets into array slice1
! ---------------------------------------
  ! Ooops - this is really expensive...
  ! basic version
  !  count=1
  !  DO cell=1,info%numtasks
  !     DO j=info%mys2,info%mye2
  !        DO k=1,info%M3
  !           DO i=info%s1(cell-1),info%e1(cell-1)
  !              out(i,j,k) = work(count) 
  !              count=count+1
  !           ENDDO
  !        ENDDO
  !     ENDDO
  !  ENDDO

  ! is this faster? Now the loops are independet from one another 
  DO cell=1,info%numtasks
     DO j=info%mys2,info%mye2
        DO k=1,info%M3
           count = 1 + (info%s1(cell-1)-1)*info%M3*(info%mye2-info%mys2+1)       &
           &      + (k-1) * (info%e1(cell-1)-info%s1(cell-1)+1)                  &
           &      + (j-info%mys2) * (info%M3*(info%e1(cell-1)-info%s1(cell-1)+1) )
           out(info%s1(cell-1):info%e1(cell-1),j,k)            &
           & = work(count:count+info%e1(cell-1)-info%s1(cell-1)) 
        ENDDO
     ENDDO
  ENDDO

! deallocate temporary memory
! ---------------------------
  DEALLOCATE(work)

END SUBROUTINE transp1d_132_1dec_to_123_2dec

