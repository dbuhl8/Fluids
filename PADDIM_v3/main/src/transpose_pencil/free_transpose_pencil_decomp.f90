!+--------------------------------------------------------------------+
!| This subroutine frees the memory allocated by the routine          |
!| init_transpose_pencil_decomposition.                               |
!+--------------------------------------------------------------------+
!| Author: Stephan Stellmach           Last modification: 29.06.06    |
!+--------------------------------------------------------------------+
SUBROUTINE free_transpose_pencil_decomp(transp_info_pencil_decomp)
  IMPLICIT NONE
  TYPE(transpose_info_pencil_decomp) :: transp_info_pencil_decomp

  CALL free_transpose_1d_decomp(transp_info_pencil_decomp%info_1st_transpose)
  CALL free_transpose_1d_decomp(transp_info_pencil_decomp%info_2nd_transpose)

END SUBROUTINE free_transpose_pencil_decomp
