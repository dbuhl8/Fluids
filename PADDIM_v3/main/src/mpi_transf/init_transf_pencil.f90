SUBROUTINE init_transf_pencil(root_communicator,nprocs1,nprocs2)
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
  INTEGER(kind=ki)         :: nprocs1,nprocs2
  INTEGER(kind=ki)         :: ierr,nprocs,myid
  LOGICAL                  :: wrong_number_of_procs
  
  ! check if nproc1*nproc2 = number of processes in root_communicator
  ! -----------------------------------------------------------------
  CALL MPI_COMM_SIZE(root_communicator,nprocs,ierr)
  CALL MPI_COMM_RANK(root_communicator,myid,ierr)
  CALL MPI_REDUCE(nprocs1*nprocs2.NE.nprocs,wrong_number_of_procs,1, &
  &               MPI_LOGICAL,MPI_LOR,0,root_communicator,ierr)
  IF ((myid.EQ.0) .AND. wrong_number_of_procs) THEN 
     PRINT*, "Error in pencil transforms initialization:"
     PRINT*, "nprocs1 * nprocs2 differs from # of processors!"
  ENDIF
  CALL MPI_BARRIER(root_communicator,ierr)

  ! initialize parallel transposes
  ! ------------------------------
  ! - since the transposes may have differing sizes in the case of zero padding,
  !   initialize two transposes, one for xy -> zx transposition, and one 
  !   for zx -> yz transposition
  CALL init_transpose_pencil_decomp(pencil_transpose_info_xy_zx,nprocs1,nprocs2, &
  &                                 Lmax+1,2*Mmax,Nz,root_communicator)
  CALL init_transpose_pencil_decomp(pencil_transpose_info_zx_yz,nprocs1,nprocs2, &
  &                                 Lmax+1,Ny,Nz,root_communicator)
  ! initialize local array bounds 
  ! -----------------------------
  ! physical space
  mysx_phys = 0
  myex_phys = Nx-1
  mysy_phys = pencil_mysy(pencil_transpose_info_zx_yz,YZ_DECOMP) - 1 
  myey_phys = pencil_myey(pencil_transpose_info_zx_yz,YZ_DECOMP) - 1 
  mysz_phys = pencil_mysz(pencil_transpose_info_zx_yz,YZ_DECOMP) - 1
  myez_phys = pencil_myez(pencil_transpose_info_zx_yz,YZ_DECOMP) - 1 

  ! spectral space
  mysx_spec = pencil_mysx(pencil_transpose_info_xy_zx,XY_DECOMP) - 1
  myex_spec = pencil_myex(pencil_transpose_info_xy_zx,XY_DECOMP) - 1
  mysy_spec = pencil_mysy(pencil_transpose_info_xy_zx,XY_DECOMP) - 1
  myey_spec = pencil_myey(pencil_transpose_info_xy_zx,XY_DECOMP) - 1
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
     WRITE(*,'(a)') '  -- y-FFTs...'
  ENDIF
  CALL init_y_transforms_pencil
  CALL MPI_BARRIER(root_communicator,ierr)
  IF (myid.EQ.0) THEN
     WRITE(*,'(a)') '     ...done'
     WRITE(*,'(a)') '  -- z-FFTs...'
  ENDIF
  CALL init_z_transforms_pencil
  CALL MPI_BARRIER(root_communicator,ierr)
  IF (myid.EQ.0) WRITE(*,'(a)') '     ...done'

END SUBROUTINE init_transf_pencil
