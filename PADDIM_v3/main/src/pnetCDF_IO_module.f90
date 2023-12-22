#include "defs_pnetCDF.h"
MODULE pnetCDF_IO_module
!+--------------------------------------------------------------------+
!| This module provides subroutines to perform the parallel IO        |
!| with parallel netCDF.                                              |
!| Currently, this is used for writing the simulation data to a disk  |
!| file for postprocession purposes (e.g. visualization) and to write |
!| restart data to a single file in a portable way. It is thus easy   |
!| to restart a computation on a different machine and with a         |
!| different number of processors.                                    |
!+--------------------------------------------------------------------+
  USE defprecision_module
  USE defs_2D_3D_module
  IMPLICIT NONE
  SAVE
!****************************************************
!****  Global Data for IO                        ****
!****************************************************

  LOGICAL  :: write_pnetCDF_sim_dat ! write simulation data in
                                    ! parallel to a netCDF file?
  LOGICAL  :: write_slice_sim_dat ! write simulation data slice to netcdf file ? 

  CHARACTER (len = 100) :: netCDF_simdat_file_name  ! file name for simulation data
  CHARACTER (len = 100) :: netCDF_xyslice_file_name  ! file name for simulation data
  CHARACTER (len = 100) :: netCDF_xzslice_file_name  ! file name for simulation data
  CHARACTER (len = 100) :: netCDF_yzslice_file_name  ! file name for simulation data

  CHARACTER (len = 100) :: netCDF_in_simdat_file_name  ! file name used in simulation data file is 
                                                       ! used to restart a simulation 
  CHARACTER (len = 100) :: netCDF_in_dump_file_name,netCDF_out_dump_file_name ! file names for restart files



  INTEGER  :: ncid_simdat            ! ncid_simdat holds the netCDF 
                                     ! id for the file containing
                                     ! (single precision) simulation data                                      

  INTEGER  :: ncid_xyslice_simdat    ! these holds the netCDF 
  INTEGER  :: ncid_xzslice_simdat    ! id for the file containingslice
  INTEGER  :: ncid_yzslice_simdat    ! (single precision) simulation data                                      

  INTEGER  :: ncid_dump              ! ncid_dump holds netCDF id for 
                                     ! restart (dump) file
  
  ! Stuff for the simulation data file 
  INTEGER  :: x_dimid_simdat,y_dimid_simdat,z_dimid_simdat, &
            & time_dimid_simdat
  INTEGER  :: x_varid_simdat,y_varid_simdat,z_varid_simdat, &
            & time_varid_simdat,timestep_varid_simdat,dt_varid_simdat
  INTEGER  :: B_therm_varid_simdat,B_comp_varid_simdat, C_Lorentz_varid_simdat ! PH
  INTEGER  :: D_visc_varid_simdat,D_therm_varid_simdat,D_comp_varid_simdat, D_mag_varid_simdat ! PH
  INTEGER  :: S_therm_varid_simdat,S_comp_varid_simdat, R_varid_simdat, Theta_varid_simdat
  INTEGER  :: Temp_varid_simdat
  INTEGER  :: Chem_varid_simdat
  INTEGER  :: ux_varid_simdat,uy_varid_simdat,uz_varid_simdat
  INTEGER  :: Bx_varid_simdat,By_varid_simdat,Bz_varid_simdat ! PH
  INTEGER  :: Gammax_varid_simdat,Gammay_varid_simdat,Gammaz_varid_simdat
  INTEGER  :: CFL_varid_simdat,dt_max_varid_simdat,dt_initial_varid_simdat
  INTEGER  :: nprocs1_varid_simdat, nprocs2_varid_simdat

  ! Stuff for the restart files 
  INTEGER  :: ri_dimid_dump,l_dimid_dump,m_dimid_dump,n_dimid_dump,time_dimid_dump
  INTEGER  :: xy_dimid_dump
  INTEGER  :: Gammax_varid_dump,Gammay_varid_dump,Gammaz_varid_dump
  INTEGER  :: B_therm_varid_dump,B_comp_varid_dump, C_Lorentz_varid_dump ! PH
  INTEGER  :: D_visc_varid_dump,D_therm_varid_dump,D_comp_varid_dump, D_mag_varid_dump ! PH
  INTEGER  :: S_therm_varid_dump,S_comp_varid_dump, R_varid_dump, Theta_varid_dump
  INTEGER  :: istep_varid_dump,time_varid_dump
  INTEGER  :: dt_varid_dump
  INTEGER  :: ri_varid_dump,kx_varid_dump,ky_varid_dump,kz_varid_dump
  INTEGER  :: Temp_varid_dump
  INTEGER  :: Chem_varid_dump
  INTEGER  :: ux_varid_dump,uy_varid_dump,uz_varid_dump
  INTEGER  :: Bx_varid_dump,By_varid_dump,Bz_varid_dump ! PH
  INTEGER  :: nprocs1_varid_dump, nprocs2_varid_dump ! DB

  ! Stuff for the XY slices
  INTEGER  :: x_dimid_xyslice,y_dimid_xyslice, &
            & time_dimid_xyslice
  INTEGER  :: x_varid_xyslice,y_varid_xyslice, &
            & time_varid_xyslice,timestep_varid_xyslice,dt_varid_xyslice
  INTEGER  :: Temp_varid_xyslice
  INTEGER  :: Chem_varid_xyslice
  INTEGER  :: ux_varid_xyslice,uy_varid_xyslice,uz_varid_xyslice
  INTEGER  :: Bx_varid_xyslice,By_varid_xyslice,Bz_varid_xyslice ! PH

  ! Stuff for the XZ slices
  INTEGER  :: x_dimid_xzslice,z_dimid_xzslice, &
            & time_dimid_xzslice
  INTEGER  :: x_varid_xzslice,z_varid_xzslice, &
            & time_varid_xzslice,timestep_varid_xzslice,dt_varid_xzslice
  INTEGER  :: Temp_varid_xzslice
  INTEGER  :: Chem_varid_xzslice
  INTEGER  :: ux_varid_xzslice,uy_varid_xzslice,uz_varid_xzslice
  INTEGER  :: Bx_varid_xzslice,By_varid_xzslice,Bz_varid_xzslice ! PH

  ! Stuff for the YZ slices
  INTEGER  :: y_dimid_yzslice,z_dimid_yzslice, &
            & time_dimid_yzslice
  INTEGER  :: y_varid_yzslice,z_varid_yzslice, &
            & time_varid_yzslice,timestep_varid_yzslice,dt_varid_yzslice
  INTEGER  :: Temp_varid_yzslice
  INTEGER  :: Chem_varid_yzslice
  INTEGER  :: ux_varid_yzslice,uy_varid_yzslice,uz_varid_yzslice
  INTEGER  :: Bx_varid_yzslice,By_varid_yzslice,Bz_varid_yzslice ! PH

CONTAINS

#include  "./pnetCDF_IO/pn_open_dump.f90"
#include  "./pnetCDF_IO/pn_write_dump.f90"
#include  "./pnetCDF_IO/pn_read_size_and_pa_from_dump.f90"
#include  "./pnetCDF_IO/pn_read_size_and_pa_from_dump_hydro.f90"
#include  "./pnetCDF_IO/pn_read_state_from_dump.f90"
#include  "./pnetCDF_IO/pn_read_state_from_dump_hydro.f90"
#include  "./pnetCDF_IO/pn_read_size_and_pa_from_simdat.f90"
#include  "./pnetCDF_IO/pn_read_state_from_simdat.f90"
#include  "./pnetCDF_IO/pn_close_dump.f90"

#include  "./pnetCDF_IO/pn_open_simdat_file.f90"
#include  "./pnetCDF_IO/pn_open_xyslice_simdat_file.f90"
#include  "./pnetCDF_IO/pn_open_xzslice_simdat_file.f90"
#include  "./pnetCDF_IO/pn_open_yzslice_simdat_file.f90"
#include  "./pnetCDF_IO/pn_write_step_simdat_file.f90"
#include  "./pnetCDF_IO/pn_write_xyslice_simdat_file.f90"
#include  "./pnetCDF_IO/pn_write_xzslice_simdat_file.f90"
#include  "./pnetCDF_IO/pn_write_yzslice_simdat_file.f90"
#include  "./pnetCDF_IO/pn_close_simdat_file.f90"
#include  "./pnetCDF_IO/pn_close_slice_simdat_file.f90"


#include  "./pnetCDF_IO/pn_check.f90"
  
END MODULE pnetCDF_IO_module
