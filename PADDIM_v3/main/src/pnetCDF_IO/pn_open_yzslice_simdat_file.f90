! This subroutine opens the data file for the simulation data stored in 
! netCDF format
SUBROUTINE pn_open_yzslice_simdat_file
  USE defprecision_module
  USE parameter_module
  USE message_passing_module, ONLY: myid
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  INTEGER :: dimids_yz(3)
  INTEGER :: date_time(8),i,j,k,ierr,myid_vert
  CHARACTER (LEN=12) :: real_clock(3)
  CHARACTER (LEN=42)  :: datestring
  REAL(kind=kr), ALLOCATABLE :: hgrid(:)
#include "pnetcdf.inc"

    IF (write_slice_sim_dat) THEN 

    ! Open the simulation data file, only for some procs.
    ! ***************************************************


    CALL pn_check( nfmpi_create(MPI_COMM_WORLD,  & 
         &         TRIM(netCDF_yzslice_file_name), &
         &         OR(NF_CLOBBER,NF_64BIT_OFFSET),MPI_INFO_NULL,ncid_yzslice_simdat                  ) )

    ! add the relevant dimensions: y,z and time 
    ! *******************************************
    CALL pn_check( nfmpi_def_dim(ncid_yzslice_simdat, "Y", INT(Ny,kind=MPI_OFFSET_KIND), &
         &         y_dimid_yzslice) )
    CALL pn_check( nfmpi_def_dim(ncid_yzslice_simdat, "Z", INT(Nz,kind=MPI_OFFSET_KIND), &
         &         z_dimid_yzslice) )
    CALL pn_check( nfmpi_def_dim(ncid_yzslice_simdat, "TIME", nfmpi_unlimited, &
         &         time_dimid_yzslice) )

    ! define variable written to netCDF file
    ! **************************************


  ! Data arrays
  
    dimids_yz = (/ y_dimid_yzslice, z_dimid_yzslice, time_dimid_yzslice /)

    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "y", NF_real,1,dimids_yz(1),y_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, y_varid_yzslice, 'long_name',  &
                 & INT(LEN('y-coordinate'),kind=MPI_OFFSET_KIND),'y-coordinate') )
    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "z", NF_real,1,dimids_yz(2),z_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, z_varid_yzslice, 'long_name',  &
                 & INT(LEN('z-coordinate'),kind=MPI_OFFSET_KIND),'z-coordinate')    )
    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "t", NF_real,1,dimids_yz(3),time_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, time_varid_yzslice, 'long_name',  &
                 & INT(LEN('time (scaled by thermal diffusion time)'),kind=MPI_OFFSET_KIND), &
                 & 'time (scaled by thermal diffusion time)'))
    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "timestep", NF_int,1,dimids_yz(3), &
                 & timestep_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, timestep_varid_yzslice,'long_name',  &
                 & INT(LEN('saved time step numbers'),kind=MPI_OFFSET_KIND),              &
                 & 'saved time step numbers'))
    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "dt", NF_real,1,dimids_yz(3), &
                 & dt_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, dt_varid_yzslice,'long_name',  &
                 & INT(LEN('time step length'),kind=MPI_OFFSET_KIND),'time step length'))



#ifdef TEMPERATURE_FIELD
    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "Temp", NF_real,3,dimids_yz,Temp_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, Temp_varid_yzslice, 'long_name',  &
                 & INT(LEN('Temperature field'),kind=MPI_OFFSET_KIND),'Temperature field')  ) 
#endif

#ifdef CHEMICAL_FIELD
    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "Chem", NF_real,3,dimids_yz,Chem_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, Chem_varid_yzslice, 'long_name',  &
                 & INT(LEN('Concentration field'),kind=MPI_OFFSET_KIND),'Concentration field') ) 
#endif

    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "ux  ", NF_real,3,dimids_yz,ux_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, ux_varid_yzslice, 'long_name',  &
                 & INT(LEN('x-component of velocity field' ),kind=MPI_OFFSET_KIND),'x-component of velocity field')) 
!#ifndef TWO_DIMENSIONAL 
    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "uy  ", NF_real,3,dimids_yz,uy_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, uy_varid_yzslice, 'long_name',  &
                 & INT(LEN('y-component of velocity field' ),kind=MPI_OFFSET_KIND),'y-component of velocity field')) 
!#endif
    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "uz  ", NF_real,3,dimids_yz,uz_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, uz_varid_yzslice, 'long_name',  &
                 & INT(LEN('z-component of velocity field' ),kind=MPI_OFFSET_KIND),'z-component of velocity field'))

#ifdef MAGNETIC
    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "Bx  ", NF_real,3,dimids_yz,Bx_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, Bx_varid_yzslice, 'long_name',  &
                 & INT(LEN('x-component of magnetic field' ),kind=MPI_OFFSET_KIND),'x-component of magnetic field')) 
!#ifndef TWO_DIMENSIONAL 
    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "By  ", NF_real,3,dimids_yz,By_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, By_varid_yzslice, 'long_name',  &
                 & INT(LEN('y-component of magnetic field' ),kind=MPI_OFFSET_KIND),'y-component of magnetic field')) 
!#endif         
    CALL pn_check( nfmpi_def_var(ncid_yzslice_simdat, "Bz  ", NF_real,3,dimids_yz,Bz_varid_yzslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, Bz_varid_yzslice, 'long_name',  &
                 & INT(LEN('z-component of magnetic field' ),kind=MPI_OFFSET_KIND),'z-component of magnetic field'))
#endif


    ! add global information 
    ! **********************
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, nf_global, 'Data',  &
                 &  INT(LEN( 'Output data of DDC code'),kind=MPI_OFFSET_KIND), &
                 &  'Output data of DDC code') )
    IF (myid.EQ.0) THEN 
       datestring = ""
       CALL DATE_AND_TIME(real_clock(1),real_clock(2),real_clock(3),date_time)
       WRITE(datestring,'(a,i4,a,i2,a,i2,a,i2,a,i2,a,i2 )')                            &
            &  "Year:",date_time(1),", month:",date_time(2),", day:",date_time(3),     &
            &  ", Time:",date_time(5),":",date_time(6),":",date_time(7)  
    ENDIF
    CALL MPI_BCAST(datestring,LEN(datestring),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL pn_check( nfmpi_put_att_text(ncid_yzslice_simdat, nf_global, 'Created',  &
                                   & INT(LEN(TRIM(datestring)),kind=MPI_OFFSET_KIND) ,TRIM(datestring)) )
    ! switch to data mode 
    ! *******************
    CALL pn_check( nfmpi_enddef(ncid_yzslice_simdat))
    ! write Parameter data and location of grid points to file 
    ! (only process 0 in independent mode) 
    CALL pn_check( nfmpi_begin_indep_data(ncid_yzslice_simdat) )
    IF (myid.EQ.0) THEN 
       ALLOCATE(hgrid(0:Ny-1))
       hgrid = (/ (j*(Gammay/REAL(Ny,kind=kr)),j=0,Ny-1) /)
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_yzslice_simdat,y_varid_yzslice,hgrid))
       DEALLOCATE(hgrid)
       ALLOCATE(hgrid(0:Nz-1))
       hgrid = (/ (k*(Gammaz/REAL(Nz,kind=kr)),k=0,Nz-1) /)
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_yzslice_simdat,z_varid_yzslice,hgrid))
       DEALLOCATE(hgrid)
       
    ENDIF
    CALL pn_check( nfmpi_end_indep_data(ncid_yzslice_simdat) )

 ENDIF
 
END SUBROUTINE pn_open_yzslice_simdat_file
