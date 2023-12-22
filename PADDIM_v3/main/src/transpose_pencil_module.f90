#include "defs_MPI.h"
!+--------------------------------------------------------------------+
!| This module provides the functionality to perform a parallel       |
!| transpose of a 3d array which is decomposed among several          |
!| processors. A 2d "pencil" decomposition is used.                   |
!+--------------------------------------------------------------------+
!| A global array of shape (N1,N2,N3) is decomposed among P processes |
!| such that each subdomain is a column of shape                      |
!|                                                                    |
!|  N1 x (N2 /~ nprocs1) x (N3 /~ nprocs2)                            |
!|                                                                    |
!| with  P = nprocs1 x nprocs2. Here, /~ denotes an "approximate"     |
!| division such that the array is decomposed among the processes     |
!| as uniformly as possible (so that every process gets columns of    |
!| roughly the sane extent in rhe x2 and x3 direction).               |
!|                                                                    |
!| The module provides everything to:                                 |
!| ----------------------------------                                 |
!|                                                                    |
!| - transpose this array from the given yz-decomposition to a        |
!|   zx-decomposition where each subdoamin is a column of shape       |
!|                                                                    |
!|           N2 x (N3 /~ nprocs2) x (N1 /~ nprocs1)                   |
!|                                                                    |
!| - and further to an xy-decomposition with subdomains of shape      |
!|                                                                    |
!|           N3 x (N1 /~ nprocs1) x (N2 /~ nprocs2).                  |
!|                                                                    |
!| - perform the above operations in reverse order                    |
!|                                                                    |
!| - compute the decomposition for given N1,N2,N3 and nprocs1,nprocs2 |
!|                                                                    |
!+--------------------------------------------------------------------+
!| Author: Stephan Stellmach                                          |
!+--------------------------------------------------------------------+
MODULE transpose_pencil_module
  USE defprecision_module
  IMPLICIT NONE
  SAVE

  PUBLIC :: init_transpose_pencil_decomp,free_transpose_pencil_decomp
  PUBLIC :: transp2d_XYZ_YZdec_to_YZX_ZXdec,transp2d_YZX_ZXdec_to_XYZ_YZdec
  PUBLIC :: transp2d_YZX_ZXdec_to_ZXY_XYdec,transp2d_ZXY_XYdec_to_YZX_ZXdec
  PUBLIC :: transp2d_XYZ_YZdec_to_XZY_ZXdec,transp2d_XZY_ZXdec_to_XYZ_YZdec
  PUBLIC :: transp2d_YZX_ZXdec_to_YXZ_XYdec,transp2d_YXZ_XYdec_to_YZX_ZXdec
  PUBLIC :: pencil_mysx,pencil_myex,pencil_mysy,pencil_myey
  PUBLIC :: pencil_mysz,pencil_myez
  PUBLIC :: transpose_info_pencil_decomp
  PUBLIC :: YZ_DECOMP,ZX_DECOMP,XY_DECOMP
  PUBLIC :: transpose_info_1d_decomp,decomp_bounds

  PUBLIC :: init_transpose_1d_decomp,free_transpose_1d_decomp
  PUBLIC :: transp1d_231_1dec_to_123_2dec,transp1d_123_2dec_to_231_1dec

  PRIVATE :: decomp_1d
  PRIVATE :: transp1d_123_2dec_to_132_1dec,transp1d_132_1dec_to_123_2dec
  
!*********************************************************************
!***                      DERIVED DATATYPES                        ***
!*********************************************************************
  TYPE transpose_info_1d_decomp
     INTEGER(kind=kiMPI)       :: communicator  ! communicator
     INTEGER(kind=ki)          :: myid          ! process id
     INTEGER(kind=ki)          :: numtasks      ! mumber of tasks
     INTEGER(kind=ki)          :: N1,N2,M3      ! The dimensions of the  
                                                ! distributed array among
                                                ! processes within communicator
                                                ! This array has the shape
                                                ! [N1,N2,M3] 
     INTEGER(kind=ki)          :: mys2,mye2     ! The local input array holds 
                                                ! global_array(1:N1,mys2:mye2,1:M3)
     INTEGER(kind=ki)          :: mys1,mye1     ! The local output array holds 
                                                ! (1:N2,1:M3,mys1:mye1)
     INTEGER(kind=ki), POINTER :: s1(:),e1(:)   ! s1(rank),e1(rank) holds mys1,mye1
                                                ! for process rank
     INTEGER(kind=ki), POINTER :: s2(:),e2(:)   ! s2(rank),e2(rank) holds mys2,mye2
                                                ! for process rank
     INTEGER(kind=ki), POINTER :: sendcounts(:) ! Arrays needed to perform
     INTEGER(kind=ki), POINTER :: recvcounts(:) ! 1d transposes with the 
     INTEGER(kind=ki), POINTER :: sdispl(:)     ! MPI_ALL_TO_ALLV
     INTEGER(kind=ki), POINTER :: rdispl(:)     ! routine
  END TYPE transpose_info_1d_decomp
  TYPE decomp_bounds
     INTEGER(kind=ki) :: mysx,myex 
     INTEGER(kind=ki) :: mysy,myey
     INTEGER(kind=ki) :: mysz,myez
  END TYPE decomp_bounds
  TYPE transpose_info_pencil_decomp
     INTEGER(kind=kiMPI) :: comm_global
     INTEGER(kind=ki)    :: nproc1,nproc2
     TYPE(decomp_bounds) :: yz_decomp
     TYPE(decomp_bounds) :: zx_decomp
     TYPE(decomp_bounds) :: xy_decomp
     TYPE(transpose_info_1d_decomp) :: info_1st_transpose
     TYPE(transpose_info_1d_decomp) :: info_2nd_transpose
  END TYPE transpose_info_pencil_decomp

!***********************
!****  Global Data  ****
!***********************

  INTEGER(kind=ki), PARAMETER :: YZ_DECOMP = 1
  INTEGER(kind=ki), PARAMETER :: ZX_DECOMP = 2
  INTEGER(kind=ki), PARAMETER :: XY_DECOMP = 3

CONTAINS

!**********************
!****  Subroutines ****
!**********************

#include "transpose_pencil/decomp_1d.f90"
#include "transpose_pencil/init_transpose_1d_decomp.f90"
#include "transpose_pencil/init_transpose_pencil_decomp.f90"
#include "transpose_pencil/free_transpose_1d_decomp.f90"
#include "transpose_pencil/free_transpose_pencil_decomp.f90"
#include "transpose_pencil/pencil_decomp_inquiry_functions.f90"
#include "transpose_pencil/create_pencil_mpi_communicators.f90"
#include "transpose_pencil/transp1d_123_2dec_to_231_1dec.f90"
#include "transpose_pencil/transp1d_231_1dec_to_123_2dec.f90"
#include "transpose_pencil/transp1d_123_2dec_to_132_1dec.f90"
#include "transpose_pencil/transp1d_132_1dec_to_123_2dec.f90"
#include "transpose_pencil/cyclic_pencil_transposes.f90"
#include "transpose_pencil/anti_cyclic_pencil_transposes.f90"

END MODULE transpose_pencil_module
