subroutine transp2d_XYZ_YZdec_to_XZY_ZXdec(xslice,yslice,transpose_info)
  IMPLICIT NONE
  COMPLEX(kind=kr)               :: xslice(:,:,:),yslice(:,:,:)
  TYPE(transpose_info_pencil_decomp) :: transpose_info

  CALL transp1d_123_2dec_to_132_1dec(xslice,yslice,                   &
  &                                  transpose_info%info_1st_transpose)

END SUBROUTINE transp2d_XYZ_YZdec_to_XZY_ZXdec


SUBROUTINE transp2d_XZY_ZXdec_to_XYZ_YZdec(yslice,xslice,transpose_info)
  IMPLICIT NONE
  COMPLEX(kind=kr)               :: xslice(:,:,:),yslice(:,:,:)
  TYPE(transpose_info_pencil_decomp) :: transpose_info

  CALL transp1d_132_1dec_to_123_2dec(yslice,xslice,                   &
  &                                  transpose_info%info_1st_transpose)

END SUBROUTINE transp2d_XZY_ZXdec_to_XYZ_YZdec


SUBROUTINE transp2d_YZX_ZXdec_to_YXZ_XYdec(yslice,zslice,transpose_info)
  IMPLICIT NONE
  COMPLEX(kind=kr)               :: yslice(:,:,:),zslice(:,:,:)
  TYPE(transpose_info_pencil_decomp) :: transpose_info

  CALL transp1d_123_2dec_to_132_1dec(yslice,zslice,                   &
  &                                  transpose_info%info_2nd_transpose)

END SUBROUTINE transp2d_YZX_ZXdec_to_YXZ_XYdec


SUBROUTINE transp2d_YXZ_XYdec_to_YZX_ZXdec(zslice,yslice,transpose_info)
  IMPLICIT NONE
  COMPLEX(kind=kr)               :: yslice(:,:,:),zslice(:,:,:)
  TYPE(transpose_info_pencil_decomp) :: transpose_info

  CALL transp1d_132_1dec_to_123_2dec(zslice,yslice,                   &
  &                                  transpose_info%info_2nd_transpose)

END SUBROUTINE transp2d_YXZ_XYdec_to_YZX_ZXdec
