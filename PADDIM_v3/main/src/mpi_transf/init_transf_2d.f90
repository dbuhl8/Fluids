SUBROUTINE init_transf_2d(root_communicator)
  USE defprecision_module
  USE transpose_pencil_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  INTEGER(kind=kiMPI)      :: root_communicator
  INTEGER(kind=kiMPI)      :: nprocs,myid
  INTEGER(kind=kiMPI)      :: ierr

 ! set nproc = number of processes in root_communicator
 ! -----------------------------------------------------------------
  CALL MPI_COMM_SIZE(root_communicator,nprocs,ierr)
  CALL MPI_COMM_RANK(root_communicator,myid,ierr)

 ! initialize parallel transposes
 ! ------------------------------
 ! - This is a bit tricky. Use the 1d transposes implemented in the transpose_pencil 
 !   module.

  CALL init_transpose_1d_decomp(two_dim_transp_info, root_communicator,Nz,Lmax+1,1)

  ! initialize local array bounds 
  ! -----------------------------
  ! physical space
  mysx_phys = 0
  myex_phys = Nx-1
  mysy_phys = 0
  myey_phys = 0
  mysz_phys = two_dim_transp_info%mys1 - 1  
  myez_phys = two_dim_transp_info%mye1 - 1  

  ! spectral space
  mysx_spec = two_dim_transp_info%mys2 - 1
  myex_spec = two_dim_transp_info%mye2 - 1  
  mysy_spec = 0
  myey_spec = 0
  mysz_spec = 0
  myez_spec = 2*Nmax-1

  ! initialize FFTW plans
  ! ---------------------
  IF (myid.EQ.0) THEN
     WRITE(*,'(a)') 'Figuring out which transform algorithm is fast'
     WRITE(*,'(a)') 'on this machine...'
     WRITE(*,'(a)') '  -- x-FFTs...'
  ENDIF
  CALL init_x_transforms_pencil
  CALL MPI_BARRIER(root_communicator,ierr)
  IF (myid.EQ.0) THEN
     WRITE(*,'(a)') '     ...done'
     WRITE(*,'(a)') '  -- z-FFTs...'
  ENDIF
  CALL init_z_transforms_pencil
  CALL MPI_BARRIER(root_communicator,ierr)
  IF (myid.EQ.0) WRITE(*,'(a)') '     ...done'

END SUBROUTINE init_transf_2d
