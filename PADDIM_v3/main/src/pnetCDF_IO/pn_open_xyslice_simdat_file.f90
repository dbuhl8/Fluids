! This subroutine opens the data file for the simulation data stored in 
! netCDF format
SUBROUTINE pn_open_xyslice_simdat_file
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
  INTEGER :: dimids_xy(3)
  INTEGER :: date_time(8),i,j,k,ierr,myid_vert
  CHARACTER (LEN=12) :: real_clock(3)
  CHARACTER (LEN=42)  :: datestring
  REAL(kind=kr), ALLOCATABLE :: hgrid(:)
#include "pnetcdf.inc"

    IF (write_slice_sim_dat) THEN 

    ! Open the simulation data file, only for some procs.
    ! ***************************************************


    CALL pn_check( nfmpi_create(MPI_COMM_WORLD,  & 
         &         TRIM(netCDF_xyslice_file_name), &
         &         OR(NF_CLOBBER,NF_64BIT_OFFSET),MPI_INFO_NULL,ncid_xyslice_simdat                  ) )

    ! add the relevant dimensions: x,y and time 
    ! *******************************************
    CALL pn_check( nfmpi_def_dim(ncid_xyslice_simdat, "X", INT(Nx,kind=MPI_OFFSET_KIND), &
         &         x_dimid_xyslice) )
    CALL pn_check( nfmpi_def_dim(ncid_xyslice_simdat, "Y", INT(Ny,kind=MPI_OFFSET_KIND), &
         &         y_dimid_xyslice) )
    CALL pn_check( nfmpi_def_dim(ncid_xyslice_simdat, "TIME", nfmpi_unlimited, &
         &         time_dimid_xyslice) )

    ! define variable written to netCDF file
    ! **************************************


  ! Data arrays
  
    dimids_xy = (/ x_dimid_xyslice, y_dimid_xyslice, time_dimid_xyslice /)

    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "x", NF_real,1,dimids_xy(1),x_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, x_varid_xyslice, 'long_name',  &
                 & INT(LEN('x-coordinate'),kind=MPI_OFFSET_KIND),'x-coordinate') )
    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "y", NF_real,1,dimids_xy(2),y_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, y_varid_xyslice, 'long_name',  &
                 & INT(LEN('y-coordinate'),kind=MPI_OFFSET_KIND),'y-coordinate')    )
    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "t", NF_real,1,dimids_xy(3),time_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, time_varid_xyslice, 'long_name',  &
                 & INT(LEN('time (scaled by thermal diffusion time)'),kind=MPI_OFFSET_KIND), &
                 & 'time (scaled by thermal diffusion time)'))
    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "timestep", NF_int,1,dimids_xy(3), &
                 & timestep_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, timestep_varid_xyslice,'long_name',  &
                 & INT(LEN('saved time step numbers'),kind=MPI_OFFSET_KIND),              &
                 & 'saved time step numbers'))
    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "dt", NF_real,1,dimids_xy(3), &
                 & dt_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, dt_varid_xyslice,'long_name',  &
                 & INT(LEN('time step length'),kind=MPI_OFFSET_KIND),'time step length'))



#ifdef TEMPERATURE_FIELD
    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "Temp", NF_real,3,dimids_xy,Temp_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, Temp_varid_xyslice, 'long_name',  &
                 & INT(LEN('Temperature field'),kind=MPI_OFFSET_KIND),'Temperature field')  ) 
#endif

#ifdef CHEMICAL_FIELD
    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "Chem", NF_real,3,dimids_xy,Chem_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, Chem_varid_xyslice, 'long_name',  &
                 & INT(LEN('Concentration field'),kind=MPI_OFFSET_KIND),'Concentration field') ) 
#endif

    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "ux  ", NF_real,3,dimids_xy,ux_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, ux_varid_xyslice, 'long_name',  &
                 & INT(LEN('x-component of velocity field' ),kind=MPI_OFFSET_KIND),'x-component of velocity field')) 
!#ifndef TWO_DIMENSIONAL 
    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "uy  ", NF_real,3,dimids_xy,uy_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, uy_varid_xyslice, 'long_name',  &
                 & INT(LEN('y-component of velocity field' ),kind=MPI_OFFSET_KIND),'y-component of velocity field')) 
!#endif
    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "uz  ", NF_real,3,dimids_xy,uz_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, uz_varid_xyslice, 'long_name',  &
                 & INT(LEN('z-component of velocity field' ),kind=MPI_OFFSET_KIND),'z-component of velocity field'))

#ifdef MAGNETIC
    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "Bx  ", NF_real,3,dimids_xy,Bx_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, Bx_varid_xyslice, 'long_name',  &
                 & INT(LEN('x-component of magnetic field' ),kind=MPI_OFFSET_KIND),'x-component of magnetic field')) 
!#ifndef TWO_DIMENSIONAL 
    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "By  ", NF_real,3,dimids_xy,By_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, By_varid_xyslice, 'long_name',  &
                 & INT(LEN('y-component of magnetic field' ),kind=MPI_OFFSET_KIND),'y-component of magnetic field')) 
!#endif         
    CALL pn_check( nfmpi_def_var(ncid_xyslice_simdat, "Bz  ", NF_real,3,dimids_xy,Bz_varid_xyslice) )
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, Bz_varid_xyslice, 'long_name',  &
                 & INT(LEN('z-component of magnetic field' ),kind=MPI_OFFSET_KIND),'z-component of magnetic field'))
#endif


    ! add global information 
    ! **********************
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, nf_global, 'Data',  &
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
    CALL pn_check( nfmpi_put_att_text(ncid_xyslice_simdat, nf_global, 'Created',  &
                                   & INT(LEN(TRIM(datestring)),kind=MPI_OFFSET_KIND) ,TRIM(datestring)) )
    ! switch to data mode 
    ! *******************
    CALL pn_check( nfmpi_enddef(ncid_xyslice_simdat))
    ! write Parameter data and location of grid points to file 
    ! (only process 0 in independent mode) 
    CALL pn_check( nfmpi_begin_indep_data(ncid_xyslice_simdat) )
    IF (myid.EQ.0) THEN 
       ALLOCATE(hgrid(0:Nx-1))
       hgrid = (/ (i*(Gammax/REAL(Nx,kind=kr)),i=0,Nx-1) /)
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_xyslice_simdat,x_varid_xyslice,hgrid))
       DEALLOCATE(hgrid)
       ALLOCATE(hgrid(0:Ny-1))
       hgrid = (/ (j*(Gammay/REAL(Ny,kind=kr)),j=0,Ny-1) /)
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_xyslice_simdat,y_varid_xyslice,hgrid))
       DEALLOCATE(hgrid)
       
    ENDIF
    CALL pn_check( nfmpi_end_indep_data(ncid_xyslice_simdat) )

 ENDIF
 
END SUBROUTINE pn_open_xyslice_simdat_file
