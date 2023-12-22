! This sumboutine adds the data of the current time step to the 
! netCDF output file 
SUBROUTINE pn_write_step_simdat_file(u,Temp,Chem,B,istep,t,dt)
  USE defprecision_module
  USE state_module
  USE mpi_transf_module, ONLY : mysy_phys,myey_phys,mysz_phys,myez_phys
  USE parameter_module
  USE message_passing_module, ONLY: myid
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  TYPE(velocity)             :: u
  TYPE(field)                :: B
  TYPE(buoyancy)             :: Temp,Chem
  REAL(kind=kr)              :: t,dt
  INTEGER(kind=ki)           :: istep
  INTEGER(kind=MPI_OFFSET_KIND)  :: starts(4),counts(4)
  INTEGER(kind=MPI_OFFSET_KIND), SAVE :: nstep_netcdf = 1
  INTEGER(kind=MPI_OFFSET_KIND), PARAMETER :: one = 1
#include "pnetcdf.inc"
  
IF (write_pnetCDF_sim_dat) THEN 
  ! write time point and time step
  ! only process 0 needs to do this...
  CALL pn_check( nfmpi_begin_indep_data(ncid_simdat) )
  IF (myid.EQ.0) THEN 
     CALL pn_check(PM_NFMPI_PUT_VARA_FLOAT(ncid_simdat, time_varid_simdat,    &
                 & nstep_netcdf, one,t             ))
     CALL pn_check(nfmpi_put_vara_int(ncid_simdat, timestep_varid_simdat,   &
                 & nstep_netcdf, one,istep        ))
     CALL pn_check(PM_NFMPI_PUT_VARA_FLOAT(ncid_simdat,dt_varid_simdat,   &
                 & nstep_netcdf, one,dt        ))
  ENDIF
  CALL pn_check( nfmpi_end_indep_data(ncid_simdat) )
  ! write physical fields in collective mode
  starts(1) = 1
  starts(2) = mysy_phys+1
  starts(3) = mysz_phys+1
  starts(4) = nstep_netcdf

  counts(1) = Nx
  counts(2) = myey_phys-mysy_phys+1
  counts(3) = myez_phys-mysz_phys+1
  counts(4) = 1

#ifdef TEMPERATURE_FIELD
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_simdat, Temp_varid_simdat,    &
                                         & starts, counts,Temp%phys)          )
#endif
#ifdef CHEMICAL_FIELD
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_simdat, Chem_varid_simdat,    &
                                         & starts, counts,Chem%phys)          )
#endif
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_simdat, ux_varid_simdat,      &
                                         & starts, counts,u%phys(:,:,:,vec_x)) )
#ifndef TWO_DIMENSIONAL
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_simdat, uy_varid_simdat,      &
                                         & starts, counts,u%phys(:,:,:,vec_y)) )
#endif
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_simdat, uz_varid_simdat,      &
                                         & starts, counts,u%phys(:,:,:,vec_z)) )

#ifdef MAGNETIC
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_simdat, Bx_varid_simdat,      &
                                         & starts, counts,B%phys(:,:,:,vec_x)) )
#ifndef TWO_DIMENSIONAL
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_simdat, By_varid_simdat,      &
                                         & starts, counts,B%phys(:,:,:,vec_y)) )
#endif
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_simdat, Bz_varid_simdat,      &
                                         & starts, counts,B%phys(:,:,:,vec_z)) )
#endif
  
  ! increase counter
  nstep_netcdf = nstep_netcdf + 1 

  ! flush all buffers 
  CALL pn_check( nfmpi_sync(ncid_simdat) )

ENDIF

END SUBROUTINE pn_write_step_simdat_file
