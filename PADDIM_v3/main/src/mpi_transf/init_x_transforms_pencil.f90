! Initialize the FFTW plans to perform r2c and c2r FFTs of the x-dependence.
! This is needed only for pencil decomposition. 
SUBROUTINE init_x_transforms_pencil
  USE defprecision_module
  USE transpose_pencil_module, ONLY : pencil_mysy,pencil_myey,  &
  &                                   pencil_mysz,pencil_myez
  IMPLICIT NONE
  INCLUDE 'fftw3.f'
  REAL(kind=kr),POINTER    :: work_real(:,:,:)
  COMPLEX(kind=kr),POINTER :: work_cmplx(:,:,:)
  INTEGER(kind=ki)         :: sy,ey,sz,ez
  INTEGER :: n,howmany
  INTEGER :: inembed,istride,idist
  INTEGER :: onembed,ostride,odist

  sy = mysy_phys 
  ey = myey_phys 
  sz = mysz_phys 
  ez = myez_phys 
  ALLOCATE(work_real(0:Nx-1,sy:ey,sz:ez))
  ALLOCATE(work_cmplx(0:Nx/2,sy:ey,sz:ez))
  ! - r2c FFT in x-direction
  n=Nx
  howmany = (ez - sz + 1) * (ey - sy + 1)
  inembed = Nx
  istride = 1
  idist = Nx
  onembed = Nx/2 + 1
  ostride = 1
  odist = Nx/2 + 1 
  CALL PM_FFTW_PLAN_MANY_DFT_R2C(plan_r2c_x_fft,1,n,howmany,         &
  &                              work_real ,inembed,istride,idist,    &
  &                              work_cmplx,onembed,ostride,odist,    &
  &                              FFTW_PATIENT+FFTW_PRESERVE_INPUT)


  ! - c2r FFT in x-direction
  n=Nx
  howmany = (ez - sz + 1) * (ey - sy + 1)
  inembed = Nx/2 + 1
  istride = 1
  idist = Nx/2 + 1 
  onembed = Nx
  ostride = 1
  odist = Nx
  CALL PM_FFTW_PLAN_MANY_DFT_C2R(plan_c2r_x_fft,1,n,howmany,         &
  &                              work_cmplx,inembed,istride,idist,   &
  &                              work_real ,onembed,ostride,odist,   &
  &                              FFTW_PATIENT+FFTW_DESTROY_INPUT)
  DEALLOCATE(work_real)
  DEALLOCATE(work_cmplx)

END SUBROUTINE init_x_transforms_pencil
