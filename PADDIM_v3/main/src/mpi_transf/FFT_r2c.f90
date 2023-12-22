! Performs a real to complex FFT for the pencil decomposition
SUBROUTINE FFT_r2c(in,out)
  USE defprecision_module
  USE transpose_pencil_module
  IMPLICIT NONE
  REAL(kind=kr)    :: in(:,:,:)
  COMPLEX(kind=kr) :: out(0:,mysx_spec:,mysy_spec:)
  COMPLEX(kind=kr), ALLOCATABLE :: xslice(:,:,:)
  COMPLEX(kind=kr), ALLOCATABLE :: work(:,:,:)
  COMPLEX(kind=kr), ALLOCATABLE :: yslice1(:,:,:),yslice2(:,:,:)
  COMPLEX(kind=kr), ALLOCATABLE :: zslice1(:,:,:),zslice2(:,:,:)
  INTEGER(kind=ki) :: i,j,k

  ! compute FFT along x 
  ALLOCATE(xslice(0:Nx/2,mysy_phys:myey_phys,mysz_phys:myez_phys))
  CALL PM_FFTW_EXECUTE_DFT_R2C(plan_r2c_x_fft,in,xslice)
#ifdef TWO_DIMENSIONAL
  ! transpose from xyz -> zxy  
  ALLOCATE(zslice1(0:Nz-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
  IF (Nx.EQ.2*Lmax) THEN 
     CALL transp1d_231_1dec_to_123_2dec(xslice,zslice1,two_dim_transp_info)
     DEALLOCATE(xslice)
  ELSE
     ALLOCATE(work(0:Lmax,mysy_phys:myey_phys,mysz_phys:myez_phys))
     ! remove modes with abs(x-wavenumber) >= Lmax 
     ! (Nyquist modes are set to zero for simplicity - they are not needed for 
     !  any calculation in the code)
     work(0:Lmax-1,:,:) = xslice(0:Lmax-1,:,:)
     work(Lmax    ,:,:) = (0._kr,0._kr)
     DEALLOCATE(xslice)
     CALL transp1d_231_1dec_to_123_2dec(work,zslice1,two_dim_transp_info)
     DEALLOCATE(work)
  ENDIF
#else
  ! transpose from xyz -> yzx
  ALLOCATE(yslice1(0:Ny-1,mysz_phys:myez_phys,mysx_spec:myex_spec))
  IF (Nx.EQ.2*Lmax) THEN 
     CALL transp2d_XYZ_YZdec_to_YZX_ZXdec(xslice,yslice1,pencil_transpose_info_zx_yz)
     DEALLOCATE(xslice)
  ELSE
     ALLOCATE(work(0:Lmax,mysy_phys:myey_phys,mysz_phys:myez_phys))
     ! remove modes with abs(x-wavenumber) >= Lmax 
     ! (Nyquist modes are set to zero for simplicity - they are not needed for 
     !  any calculation in the code)
     work(0:Lmax-1,:,:) = xslice(0:Lmax-1,:,:)
     work(Lmax    ,:,:) = (0._kr,0._kr)
     DEALLOCATE(xslice)
     CALL transp2d_XYZ_YZdec_to_YZX_ZXdec(work,yslice1,pencil_transpose_info_zx_yz)
     DEALLOCATE(work)
  ENDIF

  ! compute FFT along y 
  ALLOCATE(yslice2(0:Ny-1,mysz_phys:myez_phys,mysx_spec:myex_spec))
  CALL PM_FFTW_EXECUTE_DFT(plan_c2c_y_fft_forward,yslice1,yslice2)
  DEALLOCATE(yslice1)

  ! transpose from yzx -> zxy
  ALLOCATE(zslice1(0:Nz-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
  IF (Ny.EQ.2*Mmax) THEN
     CALL transp2d_YZX_ZXdec_to_ZXY_XYdec(yslice2,zslice1,pencil_transpose_info_xy_zx)
     DEALLOCATE(yslice2)
  ELSE
     ALLOCATE(work(0:2*Mmax-1,mysz_phys:myez_phys,mysx_spec:myex_spec))
     ! remove modes with abs(Y-wavenumber) >= mmax 
     ! (Nyquist modes are set to zero for simplicity - they are not needed for 
     !  any calculation in the code)
     DO i=mysx_spec,myex_spec
        DO k=mysz_phys,myez_phys
           work(0:Mmax-1,k,i) = yslice2(0:Mmax-1,k,i)
           work(Mmax    ,k,i) = (0._kr,0._kr)
           work(Mmax+1:2*Mmax-1,k,i) = yslice2(Ny-Mmax+1:,k,i)
        ENDDO
     ENDDO
     DEALLOCATE(yslice2)
     CALL transp2d_YZX_ZXdec_to_ZXY_XYdec(work,zslice1,pencil_transpose_info_xy_zx)
     DEALLOCATE(work)
  ENDIF
#endif

  ! compute FFT along z 
  IF (Nz.EQ.2*Nmax) THEN
     CALL PM_FFTW_EXECUTE_DFT(plan_c2c_z_fft_forward,zslice1,out)
     DEALLOCATE(zslice1)
  ELSE
     ALLOCATE(zslice2(0:Nz-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
     CALL PM_FFTW_EXECUTE_DFT(plan_c2c_z_fft_forward,zslice1,zslice2)
     DEALLOCATE(zslice1)
     ! remove modes with abs(Z-wavenumber) >= Nmax 
     ! (Nyquist modes are set to zero for simplicity - they are not needed for 
     !  any calculation in the code)
     DO j=mysy_spec,myey_spec
        DO i=mysx_spec,myex_spec
           out(0:Nmax-1,i,j) = zslice2(0:Nmax-1,i,j)
           out(Nmax    ,i,j) = (0._kr,0._kr)
           out(Nmax+1:2*Nmax-1,i,j) = zslice2(Nz-Nmax+1:,i,j)
        ENDDO
     ENDDO
     DEALLOCATE(zslice2)
  ENDIF

  ! Normalize the FFT
  out = out / (Nx*Ny*Nz)

END SUBROUTINE FFT_r2c
