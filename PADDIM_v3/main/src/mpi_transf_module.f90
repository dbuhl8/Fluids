#include "defs_MPI.h"
#include "defs_FFTW.h"
!+--------------------------------------------------------------------+
!| This module provides the functionality to perform a parallel       |
!| transform from physical space to spectral space (fourier           |
!| real-to-complex and complex to real 3D transforms).                |
!| Since the module is intended to be used for pseudospectral CFD     |
!| codes which often have to use explicit de-aliasing by the so       |
!| called 3/2-rule, the routines provide the funcitallity necessary   |
!| to implement these algorithm efficiently.                          |
!| MPI is used for explicit message passing. The FFTW 3.x library is  |
!| used to compute single processor transforms.                       |  
!|                                                                    |
!| Data layout:                                                       |
!|                                                                    |
!| In physical space, the array to be transformed is stored in the    |
!| form                                                               |
!|                                                                    |
!|       array(0:Nx-1,mysy_phys:myey_phys,mysz_phys,myez_phys),       |
!|                                                                    |
!| in spectral space                                                  |
!|                                                                    |
!|       array(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec,myey_spec).   |
!|                                                                    |       
!| The quantities mysy_phys,myey_phys,etc and mysx_spec,myex_spec,etc |
!| are computed by the routine init_transforms provided by this       |
!| module. They are public and might thus be used outside the module. |
!+--------------------------------------------------------------------+
!| subroutines:                                                       |
!|                                                                    |
!| init_transforms: initializes the transforms. Must be called before |
!|                  any transform can be computed.                    |
!| free_transforms: destroys all information collected by             |
!|                  init_transforms                                   |
!| distribute_phys: distributes data in physical space that resides   |
!|                  on just one process                               |
!| distribute_spec: distributes data in spectral space that resides   |
!|                  on just one process                               |
!| collect_phys:    collects distributed data in physical space at    |
!|                  a particular process                              |
!| collect_spec:    collects distributed data in spectral space at    |
!|                  a particular process                              |
!| FFT_r2c:         performs a 3d real to complex transform           |
!| FFT_c2r:         performs a 3d complex to real transform           |
!| collect_xy_sum_phys: Computes a sum over horizontal planes         |
!|                  in physical space.                                |
!| collect_xz_sum_phys: Computes a sum over vertical planes         |
!|                  in physical space.                                |
!+--------------------------------------------------------------------+
!| Author: Stephan Stellmach                                          |
!+--------------------------------------------------------------------+
MODULE mpi_transf_module
  USE defprecision_module
  USE transpose_pencil_module, ONLY: transpose_info_pencil_decomp,transpose_info_1d_decomp
  IMPLICIT NONE
  SAVE

  PUBLIC  :: init_transforms,FFT_r2c,FFT_c2r
  PUBLIC  :: free_transforms
  PUBLIC  :: mysx_phys,myex_phys,mysy_phys,myey_phys,mysz_phys,myez_phys
  PUBLIC  :: mysx_spec,myex_spec,mysy_spec,myey_spec,mysz_spec,myez_spec
  PUBLIC  :: collect_xy_sum_phys,collect_xz_sum_phys
  PUBLIC  :: collect_phys,collect_spec
  PUBLIC  :: distribute_phys,distribute_spec

  PRIVATE :: Nx,Ny,Nz,Lmax,Mmax,Nmax
  PRIVATE :: init_transf_pencil
  PRIVATE :: init_x_transforms_pencil,init_y_transforms_pencil,init_z_transforms_pencil
  PRIVATE :: plan_r2c_x_fft,plan_c2r_x_fft 
  PRIVATE :: plan_c2c_y_fft_forward,plan_c2c_y_fft_backward
  PRIVATE :: plan_c2c_z_fft_forward,plan_c2c_z_fft_backward
  PRIVATE :: use_FFTW_wisdom_file,FFTW_wisdom_file,uwisdom
  PRIVATE :: collect_phys_single,collect_phys_double
!
!*********************************************************************
!***                        GLOBAL DATA                            ***
!*********************************************************************

  ! local bounds of the distributed array 
  INTEGER(kind=ki) :: mysx_phys,myex_phys,mysy_phys,myey_phys,mysz_phys,myez_phys
  INTEGER(kind=ki) :: mysx_spec,myex_spec,mysy_spec,myey_spec,mysz_spec,myez_spec
  ! transform size (private to this module) 
  INTEGER(kind=ki) :: Nx,Ny,Nz,Lmax,Mmax,Nmax
  ! FFTW-Plans (private to this module)
  INTEGER*8 :: plan_r2c_x_fft,plan_c2r_x_fft 
  INTEGER*8 :: plan_c2c_y_fft_forward,plan_c2c_y_fft_backward
  INTEGER*8 :: plan_c2c_z_fft_forward,plan_c2c_z_fft_backward
  ! FFTW wisdom stuff
  LOGICAL :: use_FFTW_wisdom_file         ! Is the wisdom mechanism used?
  CHARACTER(LEN=100) :: FFTW_wisdom_file  ! Name of the FFTW wisdom file
  INTEGER(kind=ki)   :: uwisdom           ! The file unit of the wisdom file
  ! information needed to perform the transposes efficiently 
  TYPE(transpose_info_pencil_decomp)  :: pencil_transpose_info_xy_zx
  TYPE(transpose_info_pencil_decomp)  :: pencil_transpose_info_zx_yz
  TYPE(transpose_info_1d_decomp)      :: two_dim_transp_info  ! Only needed for the 2d case 
!
!*********************************************************************
!***            Interfaces for generic procedures                  ***
!*********************************************************************
INTERFACE collect_phys
   MODULE PROCEDURE collect_phys_single,collect_phys_double
END INTERFACE

CONTAINS

!*********************************************************************
!***                       Subroutines                             ***
!*********************************************************************

#include "mpi_transf/init_transforms.f90"
#include "mpi_transf/init_transf_pencil.f90"
#include "mpi_transf/init_transf_2d.f90"
#include "mpi_transf/init_x_transforms_pencil.f90"
#include "mpi_transf/init_y_transforms_pencil.f90"
#include "mpi_transf/init_z_transforms_pencil.f90"
#include "mpi_transf/FFT_r2c.f90"
#include "mpi_transf/FFT_c2r.f90"
#include "mpi_transf/free_transforms.f90"
#include "mpi_transf/f77_wisdom.f"
#include "mpi_transf/get_FFTW_wisdom.f90"
#include "mpi_transf/save_FFTW_wisdom.f90"
#include "mpi_transf/collect_phys_double.f90"
#include "mpi_transf/collect_phys_single.f90"
#include "mpi_transf/collect_spec.f90"
#include "mpi_transf/distribute_phys.f90"
#include "mpi_transf/distribute_spec.f90"
#include "mpi_transf/collect_xy_sum_phys.f90"
#include "mpi_transf/collect_xz_sum_phys.f90"
END MODULE mpi_transf_module
