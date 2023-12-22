! This subroutine reads the data from the restart file. It is assumed 
! that pn_read_size_and_pa_from_dump had already been called.
SUBROUTINE pn_read_state_from_dump_hydro(u,Temp,Chem,istep,t,dt)
  USE defprecision_module
  USE state_module, ONLY : velocity,buoyancy,ltime0
  USE parameter_module
  USE mpi_transf_module, ONLY : mysy_spec,myey_spec,mysz_spec,myez_spec,mysx_spec,myex_spec
  USE message_passing_module, ONLY : myid
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  TYPE(velocity)    :: u
  TYPE(buoyancy)    :: Temp,Chem
  INTEGER(kind=ki)  :: istep
  REAL (kind=kr)    :: t,dt
  INTEGER(kind=MPI_OFFSET_KIND)  :: starts(4),counts(4)
#include 'pnetcdf.inc'


  ! Read step, time, dt, delta = dt/dtp and deterine if this is the first step 
  CALL pn_check( nfmpi_get_var_int_all(   ncid_dump,istep_varid_dump  ,istep) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,time_varid_dump   ,t)     )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,dt_varid_dump     ,dt)    )
  PRINT*,istep,t,dt

  ! Read state from which the run will be started
  starts(1) = 1
  starts(2) = 1
  starts(3) = mysx_spec+1
  starts(4) = mysy_spec+1
  counts(1) = 2
  counts(2) = 2*Nmax
  counts(3) = myex_spec-mysx_spec+1
  counts(4) = myey_spec-mysy_spec+1

#ifdef TEMPERATURE_FIELD
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_dump, Temp_varid_dump,    &
              & starts, counts,Temp%spec(:,:,:,ltime0))                   )
#endif
#ifdef CHEMICAL_FIELD
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_dump, Chem_varid_dump,    &
              & starts, counts,Chem%spec(:,:,:,ltime0))                   )
#endif
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_dump, ux_varid_dump,    &
              & starts, counts,u%spec(:,:,:,vec_x,ltime0))              )
#ifndef TWO_DIMENSIONAL
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_dump, uy_varid_dump,    &
              & starts, counts,u%spec(:,:,:,vec_y,ltime0))              )
#endif
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_dump, uz_varid_dump,    &
              & starts, counts,u%spec(:,:,:,vec_z,ltime0))              )
  ! close dump file
  CALL pn_check( nfmpi_close(ncid_dump) )
END SUBROUTINE pn_read_state_from_dump_hydro
