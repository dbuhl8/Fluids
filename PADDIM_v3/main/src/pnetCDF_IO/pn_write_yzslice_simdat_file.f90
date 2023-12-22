!This sumboutine writes a slice of the data of the current time step to the 
! netCDF output file 
SUBROUTINE pn_write_yzslice_simdat_file(u,Temp,Chem,B,istep,t,dt)
  USE defprecision_module
  USE state_module
  USE mpi_transf_module, ONLY : mysy_phys,myey_phys,mysz_phys,myez_phys, &
        & pencil_transpose_info_zx_yz
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
  INTEGER(kind=ki)           :: istep,myid_vert,j,k
  INTEGER(kind=MPI_OFFSET_KIND)  :: starts(3),counts(3)
  INTEGER(kind=MPI_OFFSET_KIND), SAVE :: nstep_yzslice = 1
  INTEGER(kind=MPI_OFFSET_KIND), PARAMETER :: one = 1
  real(kind=kr), ALLOCATABLE :: temp_slice(:,:)
  
#include "pnetcdf.inc"
  
IF (write_slice_sim_dat) THEN 


  ! write time point and time step
  ! only process 0 needs to do this...
  CALL pn_check( nfmpi_begin_indep_data(ncid_yzslice_simdat) )
  IF (myid.EQ.0) THEN 
     CALL pn_check(PM_NFMPI_PUT_VARA_FLOAT(ncid_yzslice_simdat, time_varid_yzslice,    &
                 & nstep_yzslice, one,t             ))
     CALL pn_check(nfmpi_put_vara_int(ncid_yzslice_simdat, timestep_varid_yzslice,   &
                 & nstep_yzslice, one,istep        ))
     CALL pn_check(PM_NFMPI_PUT_VARA_FLOAT(ncid_yzslice_simdat,dt_varid_yzslice,   &
                 & nstep_yzslice, one,dt        ))
  ENDIF
  CALL pn_check( nfmpi_end_indep_data(ncid_yzslice_simdat) )
  ! write physical fields in collective mode

  ! In this case, each processor contains a little bit of the solution.
      starts(1) = mysy_phys+1
      starts(2) = mysz_phys+1
      starts(3) = nstep_yzslice
      counts(1) = myey_phys-mysy_phys+1
      counts(2) = myez_phys-mysz_phys+1
      counts(3) = 1


!  write(*,*) 'About to start', myid

   ALLOCATE(temp_slice(mysy_phys:myey_phys,mysz_phys:myez_phys))

#ifdef TEMPERATURE_FIELD
   do k=mysz_phys,myez_phys
     do j=mysy_phys,myey_phys
     temp_slice(j,k) = Temp%phys(0,j,k) 
     enddo
   enddo
!     write(*,*) 'Data copied'

   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_yzslice_simdat, Temp_varid_yzslice,    &
                                         & starts, counts,temp_slice)       )
#endif  
  

!  write(*,*) 'Printed Temp', myid
#ifdef CHEMICAL_FIELD
   do k=mysz_phys,myez_phys
     do j=mysy_phys,myey_phys
     temp_slice(j,k) = Chem%phys(0,j,k)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_yzslice_simdat, Chem_varid_yzslice,    &
                                         & starts, counts,temp_slice)          )
#endif                                 

   do k=mysz_phys,myez_phys
     do j=mysy_phys,myey_phys
     temp_slice(j,k) = u%phys(0,j,k,vec_x)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_yzslice_simdat, ux_varid_yzslice,      &
                                         & starts, counts,temp_slice) )

#ifndef TWO_DIMENSIONAL
   do k=mysz_phys,myez_phys
     do j=mysy_phys,myey_phys
     temp_slice(j,k) = u%phys(0,j,k,vec_y)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_yzslice_simdat, uy_varid_yzslice,      &
                                         & starts, counts,temp_slice) )
#endif

   do k=mysz_phys,myez_phys
     do j=mysy_phys,myey_phys
     temp_slice(j,k) = u%phys(0,j,k,vec_z)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_yzslice_simdat, uz_varid_yzslice,      &
                                         & starts, counts,temp_slice) )
 

#ifdef MAGNETIC 
   do k=mysz_phys,myez_phys
     do j=mysy_phys,myey_phys
     temp_slice(j,k) = B%phys(0,j,k,vec_x)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_yzslice_simdat, Bx_varid_yzslice,      &
                                         & starts, counts,temp_slice) )
#ifndef TWO_DIMENSIONAL
   do k=mysz_phys,myez_phys
     do j=mysy_phys,myey_phys
     temp_slice(j,k) = B%phys(0,j,k,vec_y)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_yzslice_simdat, By_varid_yzslice,      &
                                         & starts, counts,temp_slice) )
#endif
 
   do k=mysz_phys,myez_phys
     do j=mysy_phys,myey_phys
     temp_slice(j,k) = B%phys(0,j,k,vec_z)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_yzslice_simdat, Bz_varid_yzslice,      &
                                         & starts, counts,temp_slice) )
#endif
  
   DEALLOCATE(temp_slice)

!   write(*,*) ' Got here ' 

  ! increase counter
  nstep_yzslice = nstep_yzslice + 1 

  ! flush all buffers 
  CALL pn_check( nfmpi_sync(ncid_yzslice_simdat) )


ENDIF

END SUBROUTINE pn_write_yzslice_simdat_file
