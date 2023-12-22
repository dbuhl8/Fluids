SUBROUTINE pn_write_dump(u,Temp,Chem,B,t,dt,istep)
  USE defprecision_module
  USE state_module, ONLY : buoyancy,velocity,field,ltime0 ! PH
  USE message_passing_module, ONLY: myid
  USE parameter_module, ONLY: Nmax
  USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  TYPE(velocity) :: u
  TYPE(field)    :: B ! PH
  TYPE(buoyancy) :: Temp,Chem
  REAL(kind=kr) :: t,dt
  INTEGER(kind=ki) :: istep
  INTEGER(kind=MPI_OFFSET_KIND)  :: starts(4),counts(4)
#include "pnetcdf.inc"

! write current model time and time step size 
  CALL pn_check( nfmpi_begin_indep_data(ncid_dump) )
  IF (myid.EQ.0) THEN
     CALL pn_check( nfmpi_put_var_int(ncid_dump,istep_varid_dump,istep))
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,time_varid_dump,t))
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,dt_varid_dump,dt))
  ENDIF
  CALL pn_check( nfmpi_end_indep_data(ncid_dump) )

  ! write simulation data 
  starts(1) = 1
  starts(2) = 1
  starts(3) = mysx_spec+1
  starts(4) = mysy_spec+1
  counts(1) = 2
  counts(2) = 2*Nmax
  counts(3) = myex_spec-mysx_spec+1
  counts(4) = myey_spec-mysy_spec+1

#ifdef TEMPERATURE_FIELD
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_dump, Temp_varid_dump,    &
             & starts, counts,Temp%spec(:,:,:,ltime0)                )    )
#endif
#ifdef CHEMICAL_FIELD
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_dump, Chem_varid_dump,    &
             & starts, counts,Chem%spec(:,:,:,ltime0)                )    )
#endif
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_dump, ux_varid_dump,     &
             & starts, counts,u%spec(:,:,:,vec_x,ltime0)               )     )
#ifndef TWO_DIMENSIONAL
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_dump, uy_varid_dump,     &
             & starts, counts,u%spec(:,:,:,vec_y,ltime0)               )     )
#endif
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_dump, uz_varid_dump,     &
             & starts, counts,u%spec(:,:,:,vec_z,ltime0)               )     )

#ifdef MAGNETIC
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_dump, Bx_varid_dump,     &
             & starts, counts,B%spec(:,:,:,vec_x,ltime0)               )     )
#ifndef TWO_DIMENSIONAL
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_dump, By_varid_dump,     &
             & starts, counts,B%spec(:,:,:,vec_y,ltime0)               )     )
#endif
  CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_dump, Bz_varid_dump,     &
             & starts, counts,B%spec(:,:,:,vec_z,ltime0)               )     )
#endif
  
  ! flush all buffers 
  CALL pn_check( nfmpi_sync(ncid_dump) )

END SUBROUTINE pn_write_dump
