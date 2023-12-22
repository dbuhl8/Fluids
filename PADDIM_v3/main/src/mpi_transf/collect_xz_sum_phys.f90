!+---------------------------------------------------------------------+
!| This subroutine computes a vertical (xz) sum. For a global array    |
!| x(0:Nx-1,0:Ny-1,0:Nz-1) that is distributed among several processes,|
!| the quantities x(j) = sum(:,j,:) with 0 <= j < NY are computed.     |
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
SUBROUTINE collect_xz_sum_phys(local,hsum,sourcetask)
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  REAL(kind=kr)       :: local(mysx_phys:myex_phys,mysy_phys:myey_phys,mysz_phys:myez_phys)
  REAL(kind=kr)       :: hsum(0:)
  INTEGER(kind=kiMPI) :: sourcetask

  REAL(kind=kr)       :: hsum_local_xslice(mysy_phys:myey_phys)
  REAL(kind=kr)       :: hsum_local_xzplane(mysy_phys:myey_phys)
  INTEGER(kind=ki)    :: i,j,k,count,root,ierr,myid_vert,my_global_id,my_slice_id
  INTEGER(kind=kiMPI) :: sendcount
#ifdef TWO_DIMENSIONAL
  INTEGER(kind=kiMPI) :: recvcounts(two_dim_transp_info%numtasks)
  INTEGER(kind=kiMPI) :: rdispls(two_dim_transp_info%numtasks)
#else
  INTEGER(kind=kiMPI) :: recvcounts(pencil_transpose_info_zx_yz%nproc1)
  INTEGER(kind=kiMPI) :: rdispls(pencil_transpose_info_zx_yz%nproc1)
#endif

! Testing communicators
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_global_id, ierr)
  CALL MPI_COMM_RANK(pencil_transpose_info_zx_yz%info_2nd_transpose%communicator, my_slice_id, ierr)

#ifdef TWO_DIMENSIONAL
! Do nothing -- this should not be used in 2D.
#else 
  ! compute local sum in indiv. proc of all x and z points, returns y-vector
    DO j = mysy_phys,myey_phys
     hsum_local_xslice(j) = 0.
     do k = mysz_phys,myez_phys
        do i = mysx_phys,myex_phys
            hsum_local_xslice(j) = hsum_local_xslice(j) + local(i,j,k)
        enddo
     enddo
  ENDDO

  ! compute global horizontal sum in each xz-plane by summing across processors in vertical slab
  count = myey_phys - mysy_phys + 1
  root = MOD(sourcetask,pencil_transpose_info_zx_yz%nproc2) 
  CALL MPI_REDUCE(hsum_local_xslice,hsum_local_xzplane,count,                  &
  &               PM_MPI_FLOAT_TYPE,MPI_SUM,root,                           &
  &               pencil_transpose_info_zx_yz%info_2nd_transpose%communicator, &
  &               ierr)
  ! and collect everything at process sourceprocess
  myid_vert = pencil_transpose_info_zx_yz%info_2nd_transpose%myid
!  write(*,*) 'Global id is', my_global_id, 'Slice id is', my_slice_id,'myid_vert',myid_vert
  IF (myid_vert.EQ.MOD(sourcetask,pencil_transpose_info_zx_yz%nproc2)) THEN 
     sendcount = myey_phys - mysy_phys + 1 
     recvcounts(:) =  pencil_transpose_info_zx_yz%info_1st_transpose%e2(:)      &
     &              - pencil_transpose_info_zx_yz%info_1st_transpose%s2(:) + 1  
     rdispls(:) = pencil_transpose_info_zx_yz%info_1st_transpose%s2(:) - 1
     root = sourcetask / pencil_transpose_info_zx_yz%nproc2
     CALL MPI_GATHERV(hsum_local_xzplane,sendcount,PM_MPI_FLOAT_TYPE,     &
     &                hsum,recvcounts,rdispls,PM_MPI_FLOAT_TYPE,root,     &
     &                pencil_transpose_info_zx_yz%info_1st_transpose%communicator, &
     &                ierr)
  ENDIF
#endif

END SUBROUTINE collect_xz_sum_phys
