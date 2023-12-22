SUBROUTINE  allocate_fft_storage_scheme
!+--------------------------------------------------------------------+
!| This subroutine allocates memory for the matrices which contain    |
!| the FFT storage scheme.                                            |
!+--------------------------------------------------------------------+
! Allocate the matrices containing the FFT storage scheme
  IMPLICIT NONE
  ALLOCATE(kx( 0:Lmax)    )
  ALLOCATE(ky( 0:max(2*Mmax-1,0) )) 
  ALLOCATE(kz( 0:max(2*Nmax-1,0) )) 

END SUBROUTINE allocate_fft_storage_scheme
