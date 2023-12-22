SUBROUTINE  deallocate_fft_storage_scheme
!+--------------------------------------------------------------------+
!| This subroutine deallocates the arrays which are used to store     |
!| the FFT storage scheme.                                            |
!+--------------------------------------------------------------------+
  IMPLICIT NONE
  DEALLOCATE(kx,ky,kz)

END SUBROUTINE deallocate_fft_storage_scheme
