!+---------------------------------------------------------------------+
!| This subroutine computes a horizontal (xy)sum. For a global array   |
!| x(0:Nx-1,0:Ny-1,0:Nz-1) that is distributed among several processes,|
!| the quantities x(k) = sum(:,:,k) with 0 <= k < NZ are computed.     |
!| The data is assumed to be real data in physical space.              |
!|                                                                     |
!| INPUT:                                                              |
!|        local: The distributed array of values in physical space     |
!|        sourcetask: The id of the task where the sums are to be      |
!|                    collected                                        |
!|                                                                     |
!| OUTPUT:                                                             |
!|        hsum: the sum of all data in horizontal planes               |   
!|              (the result is only given on process sourcetask)       |
!|                                                                     |
!+---------------------------------------------------------------------+
SUBROUTINE collect_xy_sum_phys(local,hsum,sourcetask)
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  REAL(kind=kr)       :: local(:,:,mysz_phys:)
  REAL(kind=kr)       :: hsum(0:)
  INTEGER(kind=kiMPI) :: sourcetask

  REAL(kind=kr)       :: hsum_local_xslice(mysz_phys:myez_phys)
  REAL(kind=kr)       :: hsum_local_xyplane(mysz_phys:myez_phys)
  INTEGER(kind=ki)    :: k,count,root,ierr,myid_horiz
  INTEGER(kind=kiMPI) :: sendcount
#ifdef TWO_DIMENSIONAL
  INTEGER(kind=kiMPI) :: recvcounts(two_dim_transp_info%numtasks)
  INTEGER(kind=kiMPI) :: rdispls(two_dim_transp_info%numtasks)
#else
  INTEGER(kind=kiMPI) :: recvcounts(pencil_transpose_info_zx_yz%nproc2)
  INTEGER(kind=kiMPI) :: rdispls(pencil_transpose_info_zx_yz%nproc2)
#endif

#ifdef TWO_DIMENSIONAL
  ! compute local sum
  DO k = mysz_phys,myez_phys
     hsum_local_xyplane(k) = SUM(local(:,:,k))
  ENDDO
  sendcount = myez_phys - mysz_phys + 1
  recvcounts(:) = two_dim_transp_info%e1(:) - two_dim_transp_info%s1(:) + 1
  rdispls(:) = two_dim_transp_info%s1(:) - 1
  root = sourcetask
  CALL MPI_GATHERV(hsum_local_xyplane,sendcount,PM_MPI_FLOAT_TYPE,     &
  &                hsum,recvcounts,rdispls,PM_MPI_FLOAT_TYPE,root,     &
  &                two_dim_transp_info%communicator,ierr)
#else 
  ! compute local sum
  DO k = mysz_phys,myez_phys
     hsum_local_xslice(k) = SUM(local(:,:,k))
  ENDDO
  ! compute global horizontal sum in each xy-plane
  count = myez_phys - mysz_phys + 1
  !   - see subroutine create_pencil_MPI_communicators to understand what id the 
  !     root process has in the communicator for the 1st transpose
  root = MOD(sourcetask,pencil_transpose_info_zx_yz%nproc1) 
  CALL MPI_REDUCE(hsum_local_xslice,hsum_local_xyplane,count,                  &
  &               PM_MPI_FLOAT_TYPE,MPI_SUM,root,                           &
  &               pencil_transpose_info_zx_yz%info_1st_transpose%communicator, &
  &               ierr)
  ! and collect everything at process sourceprocess
  myid_horiz = pencil_transpose_info_zx_yz%info_1st_transpose%myid
  IF (myid_horiz.EQ. MOD(sourcetask,pencil_transpose_info_zx_yz%nproc1) ) THEN 
     sendcount = myez_phys - mysz_phys + 1 
     recvcounts(:) =  pencil_transpose_info_xy_zx%info_2nd_transpose%e2(:)      &
     &              - pencil_transpose_info_xy_zx%info_2nd_transpose%s2(:) + 1  
     rdispls(:) = pencil_transpose_info_xy_zx%info_2nd_transpose%s2(:) - 1
     root = sourcetask / pencil_transpose_info_xy_zx%nproc1
     CALL MPI_GATHERV(hsum_local_xyplane,sendcount,PM_MPI_FLOAT_TYPE,     &
     &                hsum,recvcounts,rdispls,PM_MPI_FLOAT_TYPE,root,     &
     &                pencil_transpose_info_zx_yz%info_2nd_transpose%communicator, &
     &                ierr)
  ENDIF
#endif

END SUBROUTINE collect_xy_sum_phys
