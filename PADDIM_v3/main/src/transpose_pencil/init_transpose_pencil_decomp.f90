!+--------------------------------------------------------------------+
!| This subroutine precomputes information needed to perfrom a        |
!| parallel transpose of a 3d array which is decomposed among         |
!| several processes. A 2d "pencil" decomposition is used.            |
!+--------------------------------------------------------------------+
!| A global array of shape (N1,N2,N3) is decomposed among P processes |
!| such that each subdomain is a column of size                       |
!|                                                                    |
!|  N1 x (N2 /~ nprocs1) x (N3 /~ nprocs2)                            |
!|                                                                    |
!| with  P = nprocs1 x nprocs2. Here, /~ denotes an "approximate"     |
!| division such that the array is decomposed among the processes     |
!| as uniformly as possible (so that every process gets columns of    |
!| roughly the sane extent in rhe x2 and x3 direction).               |
!| The folllowing subroutine performs several initializations which   |
!| allow to transpose this array from the given yz-decomposition      |
!| to a zx-decomposition where each subdoamin is a column of size     |
!|                                                                    |
!|  N2 x (N3 /~ nprocs2) x (N1 /~ nprocs1)                            |
!|                                                                    |
!| and further to an xy-decomposition with subdomains of size         |
!|                                                                    |
!|  N3 x (N1 /~ nprocs1) x (N2 /~ nprocs2).                           |
!|                                                                    |
!| The inverse operations also need the data initialized here.        |
!+--------------------------------------------------------------------+
!| INPUT:                                                             |
!|    nproc1,nproc2      : see above                                  |
!|    N1,N2,N3           : see above                                  |
!|    root_communicator  : the MPI communicator to be used            |
!|                         (=MPI_COMM_WORLD by default)               |
!|                                                                    |
!| OUTPUT:                                                            |
!|   transp_info_pencil_decomp : a variable of type transpose_info_   |
!|                           pencil_decomp that conains all the       |
!|                           necessary information.                   |
!+--------------------------------------------------------------------+
!| Author: Stephan Stellmach                                          |
!+--------------------------------------------------------------------+
SUBROUTINE init_transpose_pencil_decomp(transp_info_pencil_decomp,         &
&                                       nprocs1,nprocs2,N1,N2,N3,          &
&                                       root_communicator)
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  INTEGER(kind=ki)                   :: nprocs1,nprocs2,N1,N2,N3
  TYPE(transpose_info_pencil_decomp) :: transp_info_pencil_decomp
  INTEGER(kind=kiMPI),OPTIONAL       :: root_communicator
  INTEGER(kind=kiMPI)                :: basic_communicator
  INTEGER(kind=kiMPI)                :: comm_1st_transpose,comm_2nd_transpose
  INTEGER(kind=ki)                   :: myid_1st_transpose,myid_2nd_transpose
  INTEGER(kind=ki)                   :: mysx,myex,mysy1,myey1,mysy2,myey2,mysz,myez
  INTEGER(kind=ki)                   :: nproc_basic,ierr

  ! check if root_communicator is specified. If not, assume MPI_COMM_WORLD
  ! to be the root communicator 
  IF (PRESENT(root_communicator)) THEN 
     basic_communicator = root_communicator 
  ELSE
     basic_communicator = MPI_COMM_WORLD
  ENDIF
  
  transp_info_pencil_decomp%comm_global = basic_communicator

  ! check if root_communicator has nprocs1 x nprocs2 processes
  CALL MPI_COMM_SIZE(basic_communicator,nproc_basic,ierr)
  IF (nprocs1*nprocs2 .NE. nproc_basic) THEN 
     WRITE(*,'(a,I5,a)') "Root communicator does not contain ",   &
     &                   nprocs1*nprocs2," Processes."
     CALL MPI_FINALIZE(ierr)
     STOP
  ENDIF
  
  ! check if we have more processes than data points in the decomposed directions
  IF ( N1.LT.nprocs1 .OR. N2.LT.nprocs1 .OR. N2.LT.nprocs2 .OR. N3.LT.nprocs2 ) THEN
     PRINT*,"Error: To many processes for the given problem size! Aborting..."
     CALL MPI_FINALIZE(ierr)
     STOP
  ENDIF

  ! save the numbers of processes in each direction
  transp_info_pencil_decomp%nproc1 = nprocs1
  transp_info_pencil_decomp%nproc2 = nprocs2

  ! create process groups for communication 
  CALL create_pencil_MPI_communicators(nprocs1,nprocs2,                &
  &                             comm_1st_transpose,comm_2nd_transpose, &
  &                             basic_communicator)

  ! find out my id in new process groups
  CALL MPI_COMM_RANK(comm_1st_transpose,myid_1st_transpose,ierr)
  CALL MPI_COMM_RANK(comm_2nd_transpose,myid_2nd_transpose,ierr)

  ! compute domain decomposition
  CALL decomp_1d(N1,nprocs1,myid_1st_transpose,mysx,myex)
  CALL decomp_1d(N2,nprocs1,myid_1st_transpose,mysy1,myey1)
  CALL decomp_1d(N2,nprocs2,myid_2nd_transpose,mysy2,myey2)
  CALL decomp_1d(N3,nprocs2,myid_2nd_transpose,mysz,myez)

  ! save local array bounds for the decompositions 
  ! used in different stages of the program 
  !  - local array bounds for yz-decomposition
  transp_info_pencil_decomp%yz_decomp%mysx = 1
  transp_info_pencil_decomp%yz_decomp%myex = N1
  transp_info_pencil_decomp%yz_decomp%mysy = mysy1
  transp_info_pencil_decomp%yz_decomp%myey = myey1
  transp_info_pencil_decomp%yz_decomp%mysz = mysz
  transp_info_pencil_decomp%yz_decomp%myez = myez

  !  - local array bounds for zx-decomposition
  transp_info_pencil_decomp%zx_decomp%mysx = mysx
  transp_info_pencil_decomp%zx_decomp%myex = myex
  transp_info_pencil_decomp%zx_decomp%mysy = 1
  transp_info_pencil_decomp%zx_decomp%myey = N2
  transp_info_pencil_decomp%zx_decomp%mysz = mysz
  transp_info_pencil_decomp%zx_decomp%myez = myez

  !  - local array bounds for xy-decomposition
  transp_info_pencil_decomp%xy_decomp%mysx = mysx
  transp_info_pencil_decomp%xy_decomp%myex = myex
  transp_info_pencil_decomp%xy_decomp%mysy = mysy2
  transp_info_pencil_decomp%xy_decomp%myey = myey2
  transp_info_pencil_decomp%xy_decomp%mysz = 1
  transp_info_pencil_decomp%xy_decomp%myez = N3

  ! initialize transposes for 1d Decompositions
  CALL init_transpose_1d_decomp(transp_info_pencil_decomp%info_1st_transpose,   &
  &                             comm_1st_transpose,N1,N2,myez-mysz+1)
  CALL init_transpose_1d_decomp(transp_info_pencil_decomp%info_2nd_transpose,   &
  &                             comm_2nd_transpose,N2,N3,myex-mysx+1)

END SUBROUTINE init_transpose_pencil_decomp
