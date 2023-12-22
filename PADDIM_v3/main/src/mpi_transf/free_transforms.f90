!+--------------------------------------------------------------------+
!| This subroutine frees the memory where the information computed    |
!| by init_transforms is stored.                                      |
!+--------------------------------------------------------------------+
!| Author: Stephan Stellmach           Last modification: 18.09.06    |
!+--------------------------------------------------------------------+
SUBROUTINE free_transforms
  USE transpose_pencil_module
  IMPLICIT NONE

  ! delete information needed to perfrom the transposes efficiently
#ifdef TWO_DIMENSIONAL
  CALL free_transpose_1d_decomp(two_dim_transp_info)
#else
  CALL free_transpose_pencil_decomp(pencil_transpose_info_zx_yz)
  CALL free_transpose_pencil_decomp(pencil_transpose_info_xy_zx)
#endif

  ! free information computed by FFTW
  CALL PM_FFTW_CLEANUP

END SUBROUTINE free_transforms

  
