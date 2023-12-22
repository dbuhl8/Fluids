! Performs a complex to real FFT for the pencil decomposition
SUBROUTINE FFT_c2r(in,out)
  USE defprecision_module
  USE transpose_pencil_module
  IMPLICIT NONE
  COMPLEX(kind=kr) :: in(0:,mysx_spec:,mysy_spec:)
  REAL(kind=kr)    :: out(:,:,:)
  COMPLEX(kind=kr), ALLOCATABLE :: yslice1(:,:,:),yslice2(:,:,:)
  COMPLEX(kind=kr), ALLOCATABLE :: zslice(:,:,:)
  COMPLEX(kind=kr), ALLOCATABLE :: xslice(:,:,:)
  COMPLEX(kind=kr), ALLOCATABLE :: work(:,:,:)
  INTEGER(Kind=ki) :: i,j,k
  

  ! compute inverse FFT along z
  ALLOCATE(zslice(0:Nz-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
  IF (Nz.EQ.2*Nmax) THEN
     CALL PM_FFTW_EXECUTE_DFT(plan_c2c_z_fft_backward,in,zslice)
  ELSE
     ! pad z-spectrum with zeros
     ALLOCATE(work(0:Nz-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
     DO j=mysy_spec,myey_spec
        DO i=mysx_spec,myex_spec
           work(0:Nmax,i,j) = in(0:Nmax,i,j)
           work(Nmax+1:Nz-Nmax-1,i,j) = (0._kr,0._kr)
           work(Nz-Nmax:,i,j) = in(Nmax:,i,j) 
        ENDDO
     ENDDO
     CALL PM_FFTW_EXECUTE_DFT(plan_c2c_z_fft_backward,work,zslice)
     DEALLOCATE(work)
  ENDIF

#ifdef TWO_DIMENSIONAL
  ! transpose from zxy -> xyz
  ALLOCATE(xslice(0:Lmax,mysy_phys:myey_phys,mysz_phys:myez_phys))
  CALL transp1d_123_2dec_to_231_1dec(zslice,xslice,two_dim_transp_info)
#else
  ! transpose from zxy -> yzx
  ALLOCATE(yslice1(0:2*Mmax-1,mysz_phys:myez_phys,mysx_spec:myex_spec))
  CALL transp2d_ZXY_XYdec_to_YZX_ZXdec(zslice,yslice1,pencil_transpose_info_xy_zx)
  DEALLOCATE(zslice)
  ! compute inverse FFT along y 
  ALLOCATE(yslice2(0:Ny-1,mysz_phys:myez_phys,mysx_spec:myex_spec))
  IF (Ny.EQ.2*Mmax) THEN
     CALL PM_FFTW_EXECUTE_DFT(plan_c2c_y_fft_backward,yslice1,yslice2)
     DEALLOCATE(yslice1)
  ELSE
     ! pad y-spectrum with zeros
     ALLOCATE(work(0:Ny-1,mysz_phys:myez_phys,mysx_spec:myex_spec))
     DO i=mysx_spec,myex_spec
        DO k=mysz_phys,myez_phys
           work(0:Mmax,k,i) = yslice1(0:Mmax,k,i)
           work(Mmax+1:Ny-Mmax-1,k,i) = (0._kr,0._kr)
           work(Ny-Mmax:,k,i) = yslice1(Mmax:,k,i)
        ENDDO
     ENDDO
     DEALLOCATE(yslice1)
     CALL PM_FFTW_EXECUTE_DFT(plan_c2c_y_fft_backward,work,yslice2)
     DEALLOCATE(work)
  ENDIF
  ! transpose from yzx -> xyz 
  ALLOCATE(xslice(0:Lmax,mysy_phys:myey_phys,mysz_phys:myez_phys))
  CALL transp2d_YZX_ZXdec_to_XYZ_YZdec(yslice2,xslice,pencil_transpose_info_zx_yz)
  DEALLOCATE(yslice2)
#endif

  ! compute FFT along x
  IF (Nx.EQ.2*Lmax) THEN
     CALL PM_FFTW_EXECUTE_DFT_C2R(plan_c2r_x_fft,xslice,out)
     DEALLOCATE(xslice)
  ELSE
     ! pad x-spectrum with zeros
     ALLOCATE(work(0:Nx/2,mysy_phys:myey_phys,mysz_phys:myez_phys))
     DO k=mysz_phys,myez_phys
        DO j=mysy_phys,myey_phys
           work(0:Lmax,j,k) = xslice(0:Lmax,j,k)
           work(Lmax+1:,j,k) = (0._kr,0._kr)
        ENDDO
     ENDDO
     DEALLOCATE(xslice)
     CALL PM_FFTW_EXECUTE_DFT_C2R(plan_c2r_x_fft,work,out)
     DEALLOCATE(work)
  ENDIF
           
END SUBROUTINE FFT_c2r
