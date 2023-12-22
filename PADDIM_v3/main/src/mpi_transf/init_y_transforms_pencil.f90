! Initialize the FFTW plans to perform c2c FFTs of the y-dependence.
! This is needed only for pencil decomposition. 
SUBROUTINE init_y_transforms_pencil
  USE defprecision_module
  USE transpose_pencil_module, ONLY : pencil_mysx,pencil_myex,        &
  &                                          pencil_mysz,pencil_myez
  IMPLICIT NONE
  INCLUDE 'fftw3.f'
  INTEGER(kind=ki)         :: sx,ex,sz,ez
  COMPLEX(kind=kr),POINTER :: work_cmplx(:,:,:),work2_cmplx(:,:,:)
  INTEGER :: n,howmany
  INTEGER :: inembed,istride,idist
  INTEGER :: onembed,ostride,odist,sign

  sx = mysx_spec
  ex = myex_spec
  sz = mysz_phys
  ez = myez_phys
  ALLOCATE(work_cmplx(0:Ny-1,sz:ez,sx:ex))
  ALLOCATE(work2_cmplx(0:Ny-1,sz:ez,sx:ex))
  ! - forward transform
  n=Ny
  howmany = (ez - sz + 1) * (ex - sx + 1)
  inembed = Ny
  istride = 1
  idist = Ny
  onembed = Ny
  ostride = 1
  odist = Ny
  sign = FFTW_FORWARD
  CALL PM_FFTW_PLAN_MANY_DFT(plan_c2c_y_fft_forward,1,n,howmany,   &
&                            work_cmplx ,inembed,istride,idist,    &
&                            work2_cmplx,onembed,ostride,odist,    &
&                            sign,FFTW_PATIENT+FFTW_DESTROY_INPUT)
  ! - backward transform
  sign = FFTW_BACKWARD
  CALL PM_FFTW_PLAN_MANY_DFT(plan_c2c_y_fft_backward,1,n,howmany,  &
&                            work_cmplx ,inembed,istride,idist,    &
&                            work2_cmplx,onembed,ostride,odist,    &
&                            sign,FFTW_PATIENT+FFTW_DESTROY_INPUT)
  DEALLOCATE(work_cmplx)
  DEALLOCATE(work2_cmplx) 

END SUBROUTINE init_y_transforms_pencil
