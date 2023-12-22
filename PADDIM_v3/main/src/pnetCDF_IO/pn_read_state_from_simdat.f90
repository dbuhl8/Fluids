! This subroutine reads the data from the last step saved in the netcdf simulation
! data file. It is assumed that pn_read_size_and_pa_from_simdat had already been
! called.
SUBROUTINE pn_read_state_from_simdat(u,Temp,Chem,B,istep,t,dt)
  USE defprecision_module
  USE state_module, ONLY : velocity,buoyancy,field ! PH
  USE parameter_module
  USE mpi_transf_module, ONLY : mysy_phys,myey_phys,mysz_phys,myez_phys
  USE message_passing_module, ONLY : myid
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  TYPE(velocity)    :: u
  TYPE(field)       :: B ! PH
  TYPE(buoyancy)    :: Temp,Chem
  INTEGER(kind=ki)  :: istep
  INTEGER (kind=MPI_OFFSET_KIND) :: number_of_last_saved_step_simdat
  REAL (kind=kr)    :: t,dt
  INTEGER(kind=MPI_OFFSET_KIND)  :: starts(4),counts(4)
#include "pnetcdf.inc"

  CALL pn_check (nfmpi_inq_dimlen(ncid_simdat,time_dimid_simdat, &
  &              number_of_last_saved_step_simdat) )
  ! Read step, time, dt, delta = dt/dtp 
  starts(1)=number_of_last_saved_step_simdat
  counts(1)=1
  CALL pn_check( nfmpi_get_vara_int_all(   ncid_simdat,timestep_varid_simdat, &
  &                                        starts(1),counts(1),istep)         )
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_simdat,time_varid_simdat,     &
  &                                        starts(1),counts(1),t    )         )
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_simdat,dt_varid_simdat,       &
  &                                        starts(1),counts(1),dt   )         )

  ! Read state from which the run will be started
  starts(1) = 1
  starts(2) = mysy_phys+1
  starts(3) = mysz_phys+1
  starts(4) = number_of_last_saved_step_simdat

  counts(1) = Nx
  counts(2) = myey_phys-mysy_phys+1
  counts(3) = myez_phys-mysz_phys+1
  counts(4) = 1

#ifdef TEMPERATURE_FIELD
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_simdat, Temp_varid_simdat,    &
               & starts, counts,Temp%phys)                                    )
#endif
#ifdef CHEMICAL_FIELD
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_simdat, Chem_varid_simdat,    &
               & starts, counts,Chem%phys)                                    )
#endif
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_simdat, ux_varid_simdat,    &
               & starts, counts,u%phys(:,:,:,vec_x))                        )
#ifndef TWO_DIMENSIONAL
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_simdat, uy_varid_simdat,    &
               & starts, counts,u%phys(:,:,:,vec_y))                        )
#endif
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_simdat, uz_varid_simdat,    &
               & starts, counts,u%phys(:,:,:,vec_z))                        )
  
#ifdef MAGNETIC
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_simdat, Bx_varid_simdat,    &
               & starts, counts,B%phys(:,:,:,vec_x))                        )
#ifndef TWO_DIMENSIONAL
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_simdat, By_varid_simdat,    &
               & starts, counts,B%phys(:,:,:,vec_y))                        )
#endif
  CALL pn_check( PM_NFMPI_GET_VARA_FLOAT_ALL(ncid_simdat, Bz_varid_simdat,    &
               & starts, counts,B%phys(:,:,:,vec_z))                        )
#endif
  
  ! close dump file
  CALL pn_check( nfmpi_close(ncid_simdat) )
END SUBROUTINE pn_read_state_from_simdat
