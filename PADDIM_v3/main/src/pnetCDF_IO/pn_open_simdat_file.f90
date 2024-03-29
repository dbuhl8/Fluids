! This subroutine opens the data file for the simulation data stored in 
! netCDF format
SUBROUTINE pn_open_simdat_file
  USE defprecision_module
  USE parameter_module
  USE message_passing_module, ONLY: myid, nprocs1, nprocs2
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  INTEGER :: dimids(4)
  INTEGER :: date_time(8),i,j,k,ierr
  CHARACTER (LEN=12) :: real_clock(3)
  CHARACTER (LEN=42)  :: datestring
  REAL(kind=kr), ALLOCATABLE :: hgrid(:)
#include "pnetcdf.inc"

    IF (write_pnetCDF_sim_dat) THEN 

    ! Open the simulation data file
    ! *****************************

    CALL pn_check( nfmpi_create(MPI_COMM_WORLD,TRIM(netCDF_simdat_file_name), &
         &         OR(NF_CLOBBER,NF_64BIT_OFFSET),MPI_INFO_NULL,ncid_simdat                  ) )
    ! add the relevant dimensions: x,y,z and time 
    ! *******************************************
    CALL pn_check( nfmpi_def_dim(ncid_simdat, "X", INT(Nx,kind=MPI_OFFSET_KIND), &
         &         x_dimid_simdat) )
    CALL pn_check( nfmpi_def_dim(ncid_simdat, "Y", INT(Ny,kind=MPI_OFFSET_KIND), &
         &         y_dimid_simdat) )
    CALL pn_check( nfmpi_def_dim(ncid_simdat, "Z", INT(Nz,kind=MPI_OFFSET_KIND), &
         &         z_dimid_simdat) )
    CALL pn_check( nfmpi_def_dim(ncid_simdat, "TIME", nfmpi_unlimited, &
         &         time_dimid_simdat) )

    ! define variable written to netCDF file
    ! **************************************

    ! Parameter Values 
    CALL pn_check( nfmpi_def_var(ncid_simdat, "B_therm", NF_real,0,0,B_therm_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, B_therm_varid_simdat, 'long_name',  &
                 & INT(LEN('Thermal buoyancy coefficient'),kind=MPI_OFFSET_KIND),      &
                 &         'Thermal buoyancy coefficient')  )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "B_comp", NF_real,0,0,B_comp_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, B_comp_varid_simdat, 'long_name',    &
                 & INT(LEN('Compositional buoyancy coefficient'),kind=MPI_OFFSET_KIND), &
                 &         'Compositional buoyancy coefficient')  )
    ! PH{
    CALL pn_check( nfmpi_def_var(ncid_simdat, "C_Lorentz", NF_real,0,0,C_Lorentz_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, C_Lorentz_varid_simdat, 'long_name',    &
                 & INT(LEN('Lorentz force coefficient'),kind=MPI_OFFSET_KIND), &
                 &         'Lorentz force coefficient')  )
    ! }PH
    CALL pn_check( nfmpi_def_var(ncid_simdat, "D_visc", NF_real,0,0,D_visc_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, D_visc_varid_simdat, 'long_name',  &
                 & INT(LEN('Viscous diffusion parameter'),kind=MPI_OFFSET_KIND),      &
                 &         'Viscous diffusion parameter')  )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "D_therm", NF_real,0,0,D_therm_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, D_therm_varid_simdat, 'long_name',  &
                 & INT(LEN('Thermal diffusion parameter'),kind=MPI_OFFSET_KIND),       &
                 &         'Thermal diffusion parameter')  )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "D_comp", NF_real,0,0,D_comp_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, D_comp_varid_simdat, 'long_name',   &
                 & INT(LEN('Compositional diffusion parameter'),kind=MPI_OFFSET_KIND), &
                 &         'Compositional diffusion parameter')  )
    ! PH{
    CALL pn_check( nfmpi_def_var(ncid_simdat, "D_mag", NF_real,0,0,D_mag_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, D_mag_varid_simdat, 'long_name',   &
                 & INT(LEN('Magnetic diffusion parameter'),kind=MPI_OFFSET_KIND), &
                 &         'Magnetic diffusion parameter')  )
    ! }PH
    CALL pn_check( nfmpi_def_var(ncid_simdat, "S_therm", NF_real,0,0,S_therm_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, S_therm_varid_simdat, 'long_name',  &
                 & INT(LEN('Thermal background stratification parameter'),kind=MPI_OFFSET_KIND),&
                 &         'Thermal background stratification parameter')  )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "S_comp", NF_real,0,0,S_comp_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, S_comp_varid_simdat, 'long_name',    &
                 & INT(LEN('Compositional background stratification parameter'),kind=MPI_OFFSET_KIND), &
                 &         'Compositional background stratification parameter')  )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "R", NF_real,0,0,R_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, R_varid_simdat, 'long_name',    &
                 & INT(LEN('Rotational parameter'),kind=MPI_OFFSET_KIND), &
                 &         'Rotational parameter')  )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "Theta", NF_real,0,0,Theta_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, Theta_varid_simdat, 'long_name',    &
                 & INT(LEN('Angle of rotation axis and gravity parameter'),kind=MPI_OFFSET_KIND), &
                 &         'Angle of rotation axis and gravity parameter')  )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "Gammax", NF_real,0,0,Gammax_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, Gammax_varid_simdat, 'long_name',  &
                 & INT(LEN('Dimensionless length of the box'),kind=MPI_OFFSET_KIND),  &
                 & 'Dimensionless length of the box')  )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "Gammay", NF_real,0,0,Gammay_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, Gammay_varid_simdat, 'long_name',&
                 & INT(LEN('Dimensionless width of the box'),kind=MPI_OFFSET_KIND), &
                 & 'Dimensionless width of the box')  )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "Gammaz", NF_real,0,0,Gammaz_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, Gammaz_varid_simdat, 'long_name',&
                 & INT(LEN('Dimensionless height of the box'),kind=MPI_OFFSET_KIND), &
                 & 'Dimensionless height of the box')  )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "CFL", NF_real,0,0,CFL_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, CFL_varid_simdat, 'long_name',&
                 & INT(LEN('Courant-Levi-Friedrich safety factor'),kind=MPI_OFFSET_KIND),&
                 & 'Courant-Levi-Friedrich safety factor')  )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "dt_max", NF_real,0,0,dt_max_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, dt_max_varid_simdat, 'long_name',&
                 & INT(LEN('Maximum allowed time step'),kind=MPI_OFFSET_KIND),      &
                 & 'Maximum allowed time step')  )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "dt_initial", NF_real,0,0,dt_initial_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, dt_initial_varid_simdat, 'long_name',&
                 & INT(LEN('Length of first time step'),kind=MPI_OFFSET_KIND),          &
                 & 'Length of first time step')  )
    dimids = (/ x_dimid_simdat, y_dimid_simdat, z_dimid_simdat, time_dimid_simdat /)
    
    CALL pn_check( nfmpi_def_var(ncid_simdat, "x", NF_real,1,dimids(1),x_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, x_varid_simdat, 'long_name',  &
                 & INT(LEN('x-coordinate'),kind=MPI_OFFSET_KIND),'x-coordinate') )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "y", NF_real,1,dimids(2),y_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, y_varid_simdat, 'long_name',  &
                 & INT(LEN('y-coordinate'),kind=MPI_OFFSET_KIND),'y-coordinate')    )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "z", NF_real,1,dimids(3),z_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, z_varid_simdat, 'long_name',  &
                 & INT(LEN('z-coordinate'),kind=MPI_OFFSET_KIND),'z-coordinate') )
    CALL pn_check( nfmpi_def_var(ncid_simdat, "t", NF_real,1,dimids(4),time_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, time_varid_simdat, 'long_name',  &
                 & INT(LEN('time (scaled by thermal diffusion time)'),kind=MPI_OFFSET_KIND), &
                 & 'time (scaled by thermal diffusion time)'))
    CALL pn_check( nfmpi_def_var(ncid_simdat, "timestep", NF_int,1,dimids(4), &
                 & timestep_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, timestep_varid_simdat,'long_name',  &
                 & INT(LEN('saved time step numbers'),kind=MPI_OFFSET_KIND),              &
                 & 'saved time step numbers'))
    CALL pn_check( nfmpi_def_var(ncid_simdat, "dt", NF_real,1,dimids(4), &
                 & dt_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, dt_varid_simdat,'long_name',  &
                 & INT(LEN('time step length'),kind=MPI_OFFSET_KIND),'time step length'))

#ifdef TEMPERATURE_FIELD
    CALL pn_check( nfmpi_def_var(ncid_simdat, "Temp", NF_real,4,dimids,Temp_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, Temp_varid_simdat, 'long_name',  &
                 & INT(LEN('Temperature field'),kind=MPI_OFFSET_KIND),'Temperature field')  ) 
#endif
#ifdef CHEMICAL_FIELD
    CALL pn_check( nfmpi_def_var(ncid_simdat, "Chem", NF_real,4,dimids,Chem_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, Chem_varid_simdat, 'long_name',  &
                 & INT(LEN('Concentration field'),kind=MPI_OFFSET_KIND),'Concentration field') ) 
#endif
    CALL pn_check( nfmpi_def_var(ncid_simdat, "ux  ", NF_real,4,dimids,ux_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, ux_varid_simdat, 'long_name',  &
                 & INT(LEN('x-component of velocity field' ),kind=MPI_OFFSET_KIND),'x-component of velocity field')) 
    CALL pn_check( nfmpi_def_var(ncid_simdat, "uy  ", NF_real,4,dimids,uy_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, uy_varid_simdat, 'long_name',  &
                 & INT(LEN('y-component of velocity field' ),kind=MPI_OFFSET_KIND),'y-component of velocity field')) 
    CALL pn_check( nfmpi_def_var(ncid_simdat, "uz  ", NF_real,4,dimids,uz_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, uz_varid_simdat, 'long_name',  &
                 & INT(LEN('z-component of velocity field' ),kind=MPI_OFFSET_KIND),'z-component of velocity field'))

#ifdef MAGNETIC
    CALL pn_check( nfmpi_def_var(ncid_simdat, "Bx  ", NF_real,4,dimids,Bx_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, Bx_varid_simdat, 'long_name',  &
                 & INT(LEN('x-component of magnetic field' ),kind=MPI_OFFSET_KIND),'x-component of magnetic field')) 
    CALL pn_check( nfmpi_def_var(ncid_simdat, "By  ", NF_real,4,dimids,By_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, By_varid_simdat, 'long_name',  &
                 & INT(LEN('y-component of magnetic field' ),kind=MPI_OFFSET_KIND),'y-component of magnetic field')) 
    CALL pn_check( nfmpi_def_var(ncid_simdat, "Bz  ", NF_real,4,dimids,Bz_varid_simdat) )
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, Bz_varid_simdat, 'long_name',  &
                 & INT(LEN('z-component of magnetic field' ),kind=MPI_OFFSET_KIND),'z-component of magnetic field'))
#endif


    ! add global information 
    ! **********************
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, nf_global, 'Data',  &
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
    CALL pn_check( nfmpi_put_att_text(ncid_simdat, nf_global, 'Created',  &
                                   & INT(LEN(TRIM(datestring)),kind=MPI_OFFSET_KIND) ,TRIM(datestring)) )
    ! switch to data mode 
    ! *******************
    CALL pn_check( nfmpi_enddef(ncid_simdat))
    ! write Parameter data and location of grid points to file 
    ! (only process 0 in independent mode) 
    CALL pn_check( nfmpi_begin_indep_data(ncid_simdat) )
    IF (myid.EQ.0) THEN 
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,B_therm_varid_simdat,B_therm))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,B_comp_varid_simdat,B_comp))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,C_Lorentz_varid_simdat,C_Lorentz)) ! PH
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,D_visc_varid_simdat,D_visc))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,D_therm_varid_simdat,D_therm))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,D_comp_varid_simdat,D_comp))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,D_mag_varid_simdat,D_mag)) ! PH
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,S_therm_varid_simdat,S_therm))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,S_comp_varid_simdat,S_comp))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,R_varid_simdat,R))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,Theta_varid_simdat,Theta))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,Gammax_varid_simdat,Gammax))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,Gammay_varid_simdat,Gammay))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,Gammaz_varid_simdat,Gammaz))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,CFL_varid_simdat,CFL_safety_fac))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,dt_initial_varid_simdat,dt_initial))
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,dt_max_varid_simdat,dt_max))
       ALLOCATE(hgrid(0:Nx-1))
       hgrid = (/ (i*(Gammax/REAL(Nx,kind=kr)),i=0,Nx-1) /)
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,x_varid_simdat,hgrid))
       DEALLOCATE(hgrid)
       ALLOCATE(hgrid(0:Ny-1))
       hgrid = (/ (j*(Gammay/REAL(Ny,kind=kr)),j=0,Ny-1) /)
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,y_varid_simdat,hgrid))
       DEALLOCATE(hgrid)
       ALLOCATE(hgrid(0:Nz-1))
       hgrid = (/ (k*(Gammaz/REAL(Nz,kind=kr)),k=0,Nz-1) /)
       CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_simdat,z_varid_simdat,hgrid))
       DEALLOCATE(hgrid)
       
    ENDIF
    CALL pn_check( nfmpi_end_indep_data(ncid_simdat) )
 ENDIF
 
END SUBROUTINE pn_open_simdat_file
