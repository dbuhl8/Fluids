! This sumboutine writes a slice of the data of the current time step to the 
! netCDF output file 
SUBROUTINE pn_write_xzslice_simdat_file(u,Temp,Chem,B,istep,t,dt)
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
  INTEGER(kind=ki)           :: istep,myid_horiz,i,k
  INTEGER(kind=MPI_OFFSET_KIND)  :: starts(3),counts(3)
  INTEGER(kind=MPI_OFFSET_KIND), SAVE :: nstep_xzslice = 1
  INTEGER(kind=MPI_OFFSET_KIND), PARAMETER :: one = 1
  real(kind=kr), ALLOCATABLE :: temp_slice(:,:)
  
#include "pnetcdf.inc"
  
IF (write_slice_sim_dat) THEN 


  ! write time point and time step
  ! only process 0 needs to do this...
  CALL pn_check( nfmpi_begin_indep_data(ncid_xzslice_simdat) )
  IF (myid.EQ.0) THEN 
     CALL pn_check(PM_NFMPI_PUT_VARA_FLOAT(ncid_xzslice_simdat, time_varid_xzslice,    &
                 & nstep_xzslice, one,t             ))
     CALL pn_check(nfmpi_put_vara_int(ncid_xzslice_simdat, timestep_varid_xzslice,   &
                 & nstep_xzslice, one,istep        ))
     CALL pn_check(PM_NFMPI_PUT_VARA_FLOAT(ncid_xzslice_simdat,dt_varid_xzslice,   &
                 & nstep_xzslice, one,dt        ))
  ENDIF
  CALL pn_check( nfmpi_end_indep_data(ncid_xzslice_simdat) )
  ! write physical fields in collective mode

  myid_horiz = pencil_transpose_info_zx_yz%info_1st_transpose%myid
!  write(*,*) 'This proc has myid',myid,'and myid_horiz',myid_horiz

  IF(myid_horiz.eq.0) THEN  ! This should ensure that only processors with y = 0 write anything to file
      starts(1) = 1
      starts(2) = mysz_phys+1
      starts(3) = nstep_xzslice
      counts(1) = Nx
      counts(2) = myez_phys-mysz_phys+1
      counts(3) = 1
   else
      starts(1) = 1
      starts(2) = mysz_phys+1
      starts(3) = nstep_xzslice
      counts(1) = 0
      counts(2) = 0
      counts(3) = 0
   endif

!  write(*,*) 'About to start', myid

   ALLOCATE(temp_slice(0:Nx-1,mysz_phys:myez_phys))

#ifdef TEMPERATURE_FIELD
   do k=mysz_phys,myez_phys
     do i = 0,Nx-1
     temp_slice(i,k) = Temp%phys(i,mysy_phys,k) 
     enddo
   enddo
!     write(*,*) 'Data copied'

   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xzslice_simdat, Temp_varid_xzslice,    &
                                         & starts, counts,temp_slice)       )
#endif  
  

!  write(*,*) 'Printed Temp', myid
#ifdef CHEMICAL_FIELD
   do k=mysz_phys,myez_phys
     do i = 0,Nx-1
     temp_slice(i,k) = Chem%phys(i,mysy_phys,k)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xzslice_simdat, Chem_varid_xzslice,    &
                                         & starts, counts,temp_slice)          )
#endif                                 

   do k=mysz_phys,myez_phys
     do i = 0,Nx-1
     temp_slice(i,k) = u%phys(i,mysy_phys,k,vec_x)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xzslice_simdat, ux_varid_xzslice,      &
                                         & starts, counts,temp_slice) )

#ifndef TWO_DIMENSIONAL
   do k=mysz_phys,myez_phys
     do i = 0,Nx-1
     temp_slice(i,k) = u%phys(i,mysy_phys,k,vec_y)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xzslice_simdat, uy_varid_xzslice,      &
                                         & starts, counts,temp_slice) )
#endif

   do k=mysz_phys,myez_phys
     do i = 0,Nx-1
     temp_slice(i,k) = u%phys(i,mysy_phys,k,vec_z)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xzslice_simdat, uz_varid_xzslice,      &
                                         & starts, counts,temp_slice) )
  

#ifdef MAGNETIC                                 
   do k=mysz_phys,myez_phys
     do i = 0,Nx-1
     temp_slice(i,k) = B%phys(i,mysy_phys,k,vec_x)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xzslice_simdat, Bx_varid_xzslice,      &
                                         & starts, counts,temp_slice) )
#ifndef TWO_DIMENSIONAL
   do k=mysz_phys,myez_phys
     do i = 0,Nx-1
     temp_slice(i,k) = B%phys(i,mysy_phys,k,vec_y)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xzslice_simdat, By_varid_xzslice,      &
                                         & starts, counts,temp_slice) )
#endif
 
   do k=mysz_phys,myez_phys
     do i = 0,Nx-1
     temp_slice(i,k) = B%phys(i,mysy_phys,k,vec_z)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xzslice_simdat, Bz_varid_xzslice,      &
                                         & starts, counts,temp_slice) )
#endif

   DEALLOCATE(temp_slice)

!   write(*,*) ' Got here ' 

  ! increase counter
  nstep_xzslice = nstep_xzslice + 1 

  ! flush all buffers 
  CALL pn_check( nfmpi_sync(ncid_xzslice_simdat) )


ENDIF

END SUBROUTINE pn_write_xzslice_simdat_file
