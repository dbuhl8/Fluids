!+---------------------------------------------------------------------+
!| The following functions inquire which data resists on the calling   |
!| process for yz, zx and xy decomposition.                            |
!+---------------------------------------------------------------------+
!| If the global array is of the form                                  |
!|                                                                     |
!|                          a(Nx,Ny,Nz),                               |
!|                                                                     |
!| the calling process holds the data                                  |
!|                                                                     |
!| a(pencil_mysx[ , ]:pencil_myex[ , ] ,                               |
!|   pencil_mysy[ , ]:pencil_myey[ , ] ,                               |
!|   pencil_mysz[ , ]:pencil_myez[ , ])                                |
!|                                                                     |
!| The variable transpose_info, which has to be precomputed            |
!| by the routine init_transpose_pencil_decomposition, must be         |
!| provided as the first argument to pencil_mysx...pencil_myez.        |
!| The second argument is one of the three constants                   |
!|                                                                     |
!|             YZ_DECOMP, ZX_DECOMP, XY_DECOMP                         |
!|                                                                     |
!| which have been defined in the communication module and which       |
!| indicate the decomposition type (xy-, zx- or xy-decomposition)      |
!| used.                                                               |
!+---------------------------------------------------------------------+
!| Author: Stephan Stellmach            Last modification: 29.06.06    |
!+---------------------------------------------------------------------+

PURE INTEGER FUNCTION pencil_mysx(transpose_info,decomp_type)
  USE defprecision_module
  IMPLICIT NONE
  INTEGER(kind=ki), INTENT(IN) :: decomp_type
  TYPE(transpose_info_pencil_decomp), INTENT(IN) :: transpose_info

  SELECT CASE (decomp_type)
     CASE(YZ_DECOMP)
        pencil_mysx = transpose_info%yz_decomp%mysx
     CASE(ZX_DECOMP)
        pencil_mysx = transpose_info%zx_decomp%mysx
     CASE(XY_DECOMP)
        pencil_mysx = transpose_info%xy_decomp%mysx
  END SELECT
END FUNCTION pencil_mysx

PURE INTEGER FUNCTION pencil_myex(transpose_info,decomp_type)
  USE defprecision_module
  IMPLICIT NONE
  INTEGER(kind=ki), INTENT(IN) :: decomp_type
  TYPE(transpose_info_pencil_decomp), INTENT(IN) :: transpose_info

  SELECT CASE (decomp_type)
     CASE(YZ_DECOMP)
        pencil_myex = transpose_info%yz_decomp%myex
     CASE(ZX_DECOMP)
        pencil_myex = transpose_info%zx_decomp%myex
     CASE(XY_DECOMP)
        pencil_myex = transpose_info%xy_decomp%myex
  END SELECT
END FUNCTION pencil_myex

PURE INTEGER FUNCTION pencil_mysy(transpose_info,decomp_type)
  USE defprecision_module
  IMPLICIT NONE
  INTEGER(kind=ki), INTENT(IN) :: decomp_type
  TYPE(transpose_info_pencil_decomp), INTENT(IN) :: transpose_info

  SELECT CASE (decomp_type)
     CASE(YZ_DECOMP)
        pencil_mysy = transpose_info%yz_decomp%mysy
     CASE(ZX_DECOMP)
        pencil_mysy = transpose_info%zx_decomp%mysy
     CASE(XY_DECOMP)
        pencil_mysy = transpose_info%xy_decomp%mysy
  END SELECT
END FUNCTION pencil_mysy

PURE INTEGER FUNCTION pencil_myey(transpose_info,decomp_type)
  USE defprecision_module
  IMPLICIT NONE
  INTEGER(kind=ki), INTENT(IN) :: decomp_type
  TYPE(transpose_info_pencil_decomp), INTENT(IN) :: transpose_info

  SELECT CASE (decomp_type)
     CASE(YZ_DECOMP)
        pencil_myey = transpose_info%yz_decomp%myey
     CASE(ZX_DECOMP)
        pencil_myey = transpose_info%zx_decomp%myey
     CASE(XY_DECOMP)
        pencil_myey = transpose_info%xy_decomp%myey
  END SELECT
END FUNCTION pencil_myey

PURE INTEGER FUNCTION pencil_mysz(transpose_info,decomp_type)
  USE defprecision_module
  IMPLICIT NONE
  INTEGER(kind=ki), INTENT(IN) :: decomp_type
  TYPE(transpose_info_pencil_decomp), INTENT(IN) :: transpose_info

  SELECT CASE (decomp_type)
     CASE(YZ_DECOMP)
        pencil_mysz = transpose_info%yz_decomp%mysz
     CASE(ZX_DECOMP)
        pencil_mysz = transpose_info%zx_decomp%mysz
     CASE(XY_DECOMP)
        pencil_mysz = transpose_info%xy_decomp%mysz
  END SELECT
END FUNCTION pencil_mysz

PURE INTEGER FUNCTION pencil_myez(transpose_info,decomp_type)
  USE defprecision_module
  IMPLICIT NONE
  INTEGER(kind=ki), INTENT(IN) :: decomp_type
  TYPE(transpose_info_pencil_decomp), INTENT(IN) :: transpose_info

  SELECT CASE (decomp_type)
     CASE(YZ_DECOMP)
        pencil_myez = transpose_info%yz_decomp%myez
     CASE(ZX_DECOMP)
        pencil_myez = transpose_info%zx_decomp%myez
     CASE(XY_DECOMP)
        pencil_myez = transpose_info%xy_decomp%myez
  END SELECT
END FUNCTION pencil_myez






