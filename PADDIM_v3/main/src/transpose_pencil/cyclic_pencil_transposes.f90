!+--------------------------------------------------------------------+
!| The following subroutines perform parallel transposes of a 3d      |
!| complex array which is decomposed among several processes.         |
!| A 2d "pencil" decomposition is used.                               |    
!+--------------------------------------------------------------------+
!| A global array of shape (N1,N2,N3) is decomposed among P processes |
!| such that each subdomain is a column of size                       |
!|                                                                    |
!|  N1 x (N2 /~ nprocs1) x (N3 /~ nprocs2)                            |
!|                                                                    |
!| with  P = nprocs1 x nprocs2. Here, /~ denotes an "approximate"     |
!| division such that the array is decomposed among the processes     |
!| as uniformly as possible (so that every process gets columns of    |
!| roughly the sane extent in rhe x2 and x3 direction).               |
!| The folllowing subroutines transpose this array from the given     |
!| yz-decomposition to a zx-decomposition where each subdoamin is     |
!| a column of size                                                   |
!|                                                                    |
!|  N2 x (N3 /~ nprocs2) x (N1 /~ nprocs1)                            |
!|                                                                    |
!| and further to an xy-decomposition with subdomains of size         |
!|                                                                    |
!|  N3 x (N1 /~ nprocs1) x (N2 /~ nprocs2).                           |
!|                                                                    |
!| The inverse operations are also implemented.                       |
!+--------------------------------------------------------------------+
!| Author: Stephan Stellmach           Last modification: 29.06.06    |
!+--------------------------------------------------------------------+


SUBROUTINE transp2d_XYZ_YZdec_to_YZX_ZXdec(xslice,yslice,transpose_info)
  IMPLICIT NONE
  COMPLEX(kind=kr)                   :: xslice(:,:,:),yslice(:,:,:)
  TYPE(transpose_info_pencil_decomp) :: transpose_info

  CALL transp1d_123_2dec_to_231_1dec(xslice,yslice,                   &
  &                                  transpose_info%info_1st_transpose)

END SUBROUTINE transp2d_XYZ_YZdec_to_YZX_ZXdec



SUBROUTINE transp2d_YZX_ZXdec_to_XYZ_YZdec(yslice,xslice,transpose_info)
  IMPLICIT NONE
  COMPLEX(kind=kr)                   :: xslice(:,:,:),yslice(:,:,:)
  TYPE(transpose_info_pencil_decomp) :: transpose_info

  CALL transp1d_231_1dec_to_123_2dec(yslice,xslice,                   &
  &                                  transpose_info%info_1st_transpose)

END SUBROUTINE transp2d_YZX_ZXdec_to_XYZ_YZdec



SUBROUTINE transp2d_YZX_ZXdec_to_ZXY_XYdec(yslice,zslice,transpose_info)
  IMPLICIT NONE
  COMPLEX(kind=kr)                   :: yslice(:,:,:),zslice(:,:,:)
  TYPE(transpose_info_pencil_decomp) :: transpose_info

  CALL transp1d_123_2dec_to_231_1dec(yslice,zslice,                   &
  &                                  transpose_info%info_2nd_transpose)

END SUBROUTINE transp2d_YZX_ZXdec_to_ZXY_XYdec



SUBROUTINE transp2d_ZXY_XYdec_to_YZX_ZXdec(zslice,yslice,transpose_info)
  IMPLICIT NONE
  COMPLEX(kind=kr)                   :: yslice(:,:,:),zslice(:,:,:)
  TYPE(transpose_info_pencil_decomp) :: transpose_info

  CALL transp1d_231_1dec_to_123_2dec(zslice,yslice,                   &
  &                                  transpose_info%info_2nd_transpose)

END SUBROUTINE transp2d_ZXY_XYdec_to_YZX_ZXdec
