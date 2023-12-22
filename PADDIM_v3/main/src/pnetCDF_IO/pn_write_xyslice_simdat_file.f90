! This sumboutine writes a slice of the data of the current time step to the 
! netCDF output file 
SUBROUTINE pn_write_xyslice_simdat_file(u,Temp,Chem,B,istep,t,dt)
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
  INTEGER(kind=ki)           :: istep,myid_vert,i,j
  INTEGER(kind=MPI_OFFSET_KIND)  :: starts(3),counts(3)
  INTEGER(kind=MPI_OFFSET_KIND), SAVE :: nstep_xyslice = 1
  INTEGER(kind=MPI_OFFSET_KIND), PARAMETER :: one = 1
  real(kind=kr), ALLOCATABLE :: temp_slice(:,:)
  
#include "pnetcdf.inc"
  
IF (write_slice_sim_dat) THEN 


  ! write time point and time step
  ! only process 0 needs to do this...
  CALL pn_check( nfmpi_begin_indep_data(ncid_xyslice_simdat) )
  IF (myid.EQ.0) THEN 
     CALL pn_check(PM_NFMPI_PUT_VARA_FLOAT(ncid_xyslice_simdat, time_varid_xyslice,    &
                 & nstep_xyslice, one,t             ))
     CALL pn_check(nfmpi_put_vara_int(ncid_xyslice_simdat, timestep_varid_xyslice,   &
                 & nstep_xyslice, one,istep        ))
     CALL pn_check(PM_NFMPI_PUT_VARA_FLOAT(ncid_xyslice_simdat,dt_varid_xyslice,   &
                 & nstep_xyslice, one,dt        ))
  ENDIF
  CALL pn_check( nfmpi_end_indep_data(ncid_xyslice_simdat) )
  ! write physical fields in collective mode

  myid_vert = pencil_transpose_info_zx_yz%info_2nd_transpose%myid
!  write(*,*) 'This proc has myid',myid,'and myid_vert',myid_vert

  IF(myid_vert.eq.0) THEN  ! This should ensure that only processors with z = 0 write anything to file
      starts(1) = 1
      starts(2) = mysy_phys+1
      starts(3) = nstep_xyslice
      counts(1) = Nx
      counts(2) = myey_phys-mysy_phys+1
      counts(3) = 1

   else
      starts(1) = 1
      starts(2) = mysy_phys+1
      starts(3) = nstep_xyslice
      counts(1) = 0
      counts(2) = 0
      counts(3) = 0
   endif

!  write(*,*) 'About to start', myid

   ALLOCATE(temp_slice(0:Nx-1,mysy_phys:myey_phys))

#ifdef TEMPERATURE_FIELD
   do j=mysy_phys,myey_phys
     do i = 0,Nx-1
     temp_slice(i,j) = Temp%phys(i,j,mysz_phys) 
     enddo
   enddo
!     write(*,*) 'Data copied'

   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xyslice_simdat, Temp_varid_xyslice,    &
                                         & starts, counts,temp_slice)       )
#endif  
  

!  write(*,*) 'Printed Temp', myid
#ifdef CHEMICAL_FIELD
   do j=mysy_phys,myey_phys
     do i = 0,Nx-1
     temp_slice(i,j) = Chem%phys(i,j,mysz_phys)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xyslice_simdat, Chem_varid_xyslice,    &
                                         & starts, counts,temp_slice)          )
#endif                                 

   do j=mysy_phys,myey_phys
     do i = 0,Nx-1
     temp_slice(i,j) = u%phys(i,j,mysz_phys,vec_x)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xyslice_simdat, ux_varid_xyslice,      &
                                         & starts, counts,temp_slice) )

#ifndef TWO_DIMENSIONAL
   do j=mysy_phys,myey_phys
     do i = 0,Nx-1
     temp_slice(i,j) = u%phys(i,j,mysz_phys,vec_y)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xyslice_simdat, uy_varid_xyslice,      &
                                         & starts, counts,temp_slice) )
#endif

   do j=mysy_phys,myey_phys
     do i = 0,Nx-1
     temp_slice(i,j) = u%phys(i,j,mysz_phys,vec_z)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xyslice_simdat, uz_varid_xyslice,      &
                                         & starts, counts,temp_slice) )

#ifdef MAGNETIC
   do j=mysy_phys,myey_phys
     do i = 0,Nx-1
     temp_slice(i,j) = B%phys(i,j,mysz_phys,vec_x)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xyslice_simdat, Bx_varid_xyslice,      &
                                         & starts, counts,temp_slice) )
#ifndef TWO_DIMENSIONAL
   do j=mysy_phys,myey_phys
     do i = 0,Nx-1
     temp_slice(i,j) = B%phys(i,j,mysz_phys,vec_y)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xyslice_simdat, By_varid_xyslice,      &
                                         & starts, counts,temp_slice) )
#endif
 
   do j=mysy_phys,myey_phys
     do i = 0,Nx-1
     temp_slice(i,j) = B%phys(i,j,mysz_phys,vec_z)
     enddo
   enddo
   CALL pn_check( PM_NFMPI_PUT_VARA_FLOAT_ALL(ncid_xyslice_simdat, Bz_varid_xyslice,      &
                                         & starts, counts,temp_slice) )

#endif
   DEALLOCATE(temp_slice)

!   write(*,*) ' Got here ' 

  ! increase counter
  nstep_xyslice = nstep_xyslice + 1 

  ! flush all buffers 
  CALL pn_check( nfmpi_sync(ncid_xyslice_simdat) )


ENDIF

END SUBROUTINE pn_write_xyslice_simdat_file
