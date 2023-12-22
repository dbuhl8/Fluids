! The follwing subroutine transposes a three dimensional 
! array in(:,:,:) which is decomposed along the x2-direction to an 
! array out(:,:,:) which is decomposed along the x1-direction. 
! The order in which the different directions are stored changes 
! from x1:x2:x3 -> x1:x3:x2.  
SUBROUTINE transp1d_123_2dec_to_132_1dec(in,out,info)
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  TYPE(transpose_info_1d_decomp) :: info
  COMPLEX(kind=kr)               :: in(:,info%mys2:,:)
  COMPLEX(kind=kr)               :: out(info%mys1:,:,:)
  COMPLEX(kind=kr),ALLOCATABLE   :: work(:)
  INTEGER(kind=ki)               :: ierr
  INTEGER(kind=ki)               :: j,k,cell
  INTEGER(kind=kli)              :: count,work_length
! allocate temporary memory
! ------------------------- 
  work_length =  INT(info%mye2-info%mys2+1,kind=kli) &
  &             *INT(info%N1,kind=kli)               &
  &             *INT(info%M3,kind=kli)
  ALLOCATE(work(work_length))
! sort send packets into continous parts of array work
! ----------------------------------------------------
  ! Ooops - this is really expensive...
  ! basic version
  !  count=1
  !  DO cell=1,info%numtasks
  !     DO j=info%mys2,info%mye2
  !        DO k=1,info%M3
  !           DO i=info%s1(cell-1),info%e1(cell-1)
  !              work(count) = in(i,j,k)
  !              count=count+1
  !           ENDDO
  !        ENDDO
  !     ENDDO
  !  ENDDO

  ! is this faster? Now the loops are independet
  DO cell=1,info%numtasks
     DO j=info%mys2,info%mye2
        DO k=1,info%M3
           count = 1 + (info%s1(cell-1)-1)*info%M3*(info%mye2-info%mys2+1)       &
           &      + (k-1) * (info%e1(cell-1)-info%s1(cell-1)+1)                  &
           &      + (j-info%mys2) * (info%M3*(info%e1(cell-1)-info%s1(cell-1)+1) )
           work(count:count+info%e1(cell-1)-info%s1(cell-1)) &
           & = in(info%s1(cell-1):info%e1(cell-1),j,k)
        ENDDO
     ENDDO
  ENDDO

  ! or this?
  !FORALL(cell=1:info%numtasks,j=info%mys2:info%mye2,k=1:info%M3)            &
  !& work(  1 + (info%s1(cell-1)-1)*info%M3*(info%mye2-info%mys2+1)          &
  !&        + (k-1) * (info%e1(cell-1)-info%s1(cell-1)+1)                    &
  !&        + (j-info%mys2) * (info%M3*(info%e1(cell-1)-info%s1(cell-1)+1))  &
  !&      : 1 + (info%s1(cell-1)-1)*info%M3*(info%mye2-info%mys2+1)          &
  !&        + (k-1) * (info%e1(cell-1)-info%s1(cell-1)+1)                    &
  !&        + (j-info%mys2) * (info%M3*(info%e1(cell-1)-info%s1(cell-1)+1))  &
  !&        +info%e1(cell-1)-info%s1(cell-1)                               ) &
  !& = in(info%s1(cell-1):info%e1(cell-1),j,k)

     

! perform the all to all communication using MPI_ALLTOALLV
! --------------------------------------------------------  
  CALL MPI_ALLTOALLV(work,info%sendcounts,info%sdispl,PM_MPI_COMPLEX_TYPE,    &
  &                  out ,info%recvcounts,info%rdispl,PM_MPI_COMPLEX_TYPE,    &
  &                  info%communicator,ierr)
! deallocate temporary memory 
! ---------------------------
  DEALLOCATE(work)

END SUBROUTINE transp1d_123_2dec_to_132_1dec
