! Initialize the FFTW plans to perform c2c FFTs of the z-dependence.
! This is needed only for pencil decomposition. 
SUBROUTINE init_z_transforms_pencil
  USE defprecision_module
  USE transpose_pencil_module, ONLY : pencil_mysx,pencil_myex, &
  &                                   pencil_mysz,pencil_myez
  IMPLICIT NONE
  INCLUDE 'fftw3.f'
  INTEGER(kind=ki)         :: sx,ex,sy,ey
  COMPLEX(kind=kr),POINTER :: work_cmplx(:,:,:),work2_cmplx(:,:,:)
  INTEGER :: n,howmany
  INTEGER :: inembed,istride,idist
  INTEGER :: onembed,ostride,odist,sign

  sx = mysx_spec
  ex = myex_spec
  sy = mysy_spec
  ey = myey_spec
  ALLOCATE(work_cmplx(0:Nz-1,sx:ex,sy:ey))
  ALLOCATE(work2_cmplx(0:Nz-1,sx:ex,sy:ey))
  ! - forward transform
  n=Nz
  howmany = (ey - sy + 1) * (ex - sx + 1)
  inembed = Nz
  istride = 1
  idist = Nz
  onembed = Nz
  ostride = 1
  odist = Nz
  sign = FFTW_FORWARD
  CALL PM_FFTW_PLAN_MANY_DFT(plan_c2c_z_fft_forward,1,n,howmany,   &
&                            work_cmplx ,inembed,istride,idist,    &
&                            work2_cmplx,onembed,ostride,odist,    &
&                            sign,FFTW_PATIENT+FFTW_DESTROY_INPUT)
  ! - backward transform
  sign = FFTW_BACKWARD
  CALL PM_FFTW_PLAN_MANY_DFT(plan_c2c_z_fft_backward,1,n,howmany,  &
&                            work_cmplx ,inembed,istride,idist,    &
&                            work2_cmplx,onembed,ostride,odist,    &
&                            sign,FFTW_PATIENT+FFTW_PRESERVE_INPUT)
  DEALLOCATE(work_cmplx)
  DEALLOCATE(work2_cmplx) 

END SUBROUTINE init_z_transforms_pencil
