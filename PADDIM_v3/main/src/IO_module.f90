#include "defs_MPI.h"
#include "../../stuff_needed/include/jcmagic.h"
MODULE IO_module
!+--------------------------------------------------------------------+
!| This module provides subroutines to perform the necessary file IO. |
!+--------------------------------------------------------------------+
 USE defprecision_module
 USE defs_2D_3D_module
 USE parameter_module, ONLY:
 USE message_passing_module, ONLY:
 USE state_module, ONLY:
 USE testing_module, ONLY:
 USE pnetCDF_IO_module, ONLY:
 IMPLICIT NONE
 SAVE
 !****************************************************
 !****  Global Data for IO                        ****
 !****************************************************
 INTEGER (kind=ki),PARAMETER :: nout=6              ! Number of formatted output files
 INTEGER (kind=ki)           :: uout(nout)          ! Array containing the output file units
 CHARACTER (LEN=100)         :: outfile(nout)       ! Array containing the output file names 
 CHARACTER (LEN=100)         :: jc_out              ! File that contains compressed fields
 INTEGER (kind=ki)           :: number_of_jc_fields ! Number of scalar fields that are to be
                                                    ! written to JC-Files
 LOGICAL                     :: write_compressed_fields   ! should physical fields be written
                                                          ! to jpeg compressed output files?
 INTEGER(kind=ki) :: n_comp_diag  ! Compute diagnostics every n-th time step
 INTEGER(kind=ki) :: n_wrt_jc     ! Write jc-file every n-th time step 
 INTEGER(kind=ki) :: n_wrt_netCDF ! Write netCFD file every n-th time step
 INTEGER(kind=ki) :: n_wrt_slice  ! Write slice file every n-th time step
 INTEGER(kind=ki) :: n_wrt_dump   ! Write restart file every n-th time step 
 INTEGER(kind=ki) :: n_wrt_spec   ! Write spectra file every n-th time step 
 INTEGER(kind=ki) :: n_wrt_prof   ! Write profile file every n-th time step 
CONTAINS

#include "IO/open_files.f90"
#include "IO/write_output_files.f90"
#include "IO/write_compressed_file.f90"
#include "IO/write_diagnostics_file.f90"
#include "IO/write_horizontal_spectra.f90"
#include "IO/write_vertical_spectra.f90"
#include "IO/write_shell_spectra.f90"
#include "IO/write_z_profile.f90"
#include "IO/write_y_profile.f90"
#include "IO/close_files.f90"

END MODULE IO_module
