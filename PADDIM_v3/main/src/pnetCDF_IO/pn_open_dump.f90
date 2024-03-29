SUBROUTINE pn_open_dump
  USE defprecision_module
  USE parameter_module, ONLY : Nmax,Lmax,Mmax,B_therm,B_comp,C_Lorentz,D_visc,D_therm,D_comp,D_mag, &
  &                            S_therm,S_comp,R,Theta,Gammax,Gammay,Gammaz,kx,ky,kz
  USE message_passing_module, ONLY: myid, nprocs1, nprocs2
#ifdef MPI_MODULE
  USE MPI 
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
#include "pnetcdf.inc" 
  INTEGER :: dimids(4),ierr,i
  CHARACTER (LEN=12) :: real_clock(3)
  INTEGER           :: date_time(8)
  CHARACTER (LEN=42)  :: datestring
  
  ! Open the restart file - if it exists aready, just overwrite it
  ! **************************************************************
  CALL pn_check( nfmpi_create(MPI_COMM_WORLD,TRIM(netCDF_out_dump_file_name), &
       &         OR(NF_CLOBBER,NF_64BIT_OFFSET),MPI_INFO_NULL,ncid_dump                  ) )
 
  ! add the relevant dimensions
  ! ***************************
  CALL pn_check( nfmpi_def_dim(ncid_dump, "ri", INT(2,  kind=MPI_OFFSET_KIND), &
       &         ri_dimid_dump) )
  CALL pn_check( nfmpi_def_dim(ncid_dump, "n", INT(2*Nmax,  kind=MPI_OFFSET_KIND), &
       &         n_dimid_dump) )
  CALL pn_check( nfmpi_def_dim(ncid_dump, "l", INT(Lmax+1,  kind=MPI_OFFSET_KIND), &
       &         l_dimid_dump) )
  CALL pn_check( nfmpi_def_dim(ncid_dump, "m", INT(MAX(2*Mmax,1),  kind=MPI_OFFSET_KIND), &
       &         m_dimid_dump) ) !in 2D, the value m=0 is stored as well.

  ! define variables written to netCDF file
  ! ***************************************

  ! Parameter Values
  CALL pn_check( nfmpi_def_var(ncid_dump, "Gammax", PM_NF_FLOAT,0,0,Gammax_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, Gammax_varid_dump, 'long_name',&
               & INT(LEN('Dimensionless length of the box'),kind=MPI_OFFSET_KIND), &
               & 'Dimensionless length of the box')  )
#ifndef TWO_DIMENSIONAL
  CALL pn_check( nfmpi_def_var(ncid_dump, "Gammay", PM_NF_FLOAT,0,0,Gammay_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, Gammay_varid_dump, 'long_name',&
               & INT(LEN('Dimensionless width of the box'),kind=MPI_OFFSET_KIND),  &
               & 'Dimensionless width of the box')  )
#endif
  CALL pn_check( nfmpi_def_var(ncid_dump, "Gammaz", PM_NF_FLOAT,0,0,Gammaz_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, Gammaz_varid_dump, 'long_name',&
               & INT(LEN('Dimensionless height of the box'),kind=MPI_OFFSET_KIND),  &
               & 'Dimensionless height of the box')  )
  CALL pn_check( nfmpi_def_var(ncid_dump, "B_therm", PM_NF_FLOAT,0,0,B_therm_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, B_therm_varid_dump, 'long_name',  &
               & INT(LEN('Thermal buoyancy coefficient'),kind=MPI_OFFSET_KIND), &
               & 'Thermal buoyancy coefficient')  )
  CALL pn_check( nfmpi_def_var(ncid_dump, "B_comp", PM_NF_FLOAT,0,0,B_comp_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, B_comp_varid_dump, 'long_name',  &
               & INT(LEN('Compositional buoyancy coefficient'),kind=MPI_OFFSET_KIND),  &
               & 'Compositional buoyancy coefficient')  )
  ! PH{
  CALL pn_check( nfmpi_def_var(ncid_dump, "C_Lorentz", PM_NF_FLOAT,0,0,C_Lorentz_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, C_Lorentz_varid_dump, 'long_name',  &
               & INT(LEN('Lorentz force coefficient'),kind=MPI_OFFSET_KIND),  &
               & 'Lorentz force coefficient')  )
  ! }PH
  CALL pn_check( nfmpi_def_var(ncid_dump, "D_visc", PM_NF_FLOAT,0,0,D_visc_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump,D_visc_varid_dump, 'long_name',  &
               & INT(LEN('Viscous diffusion parameter'),kind=MPI_OFFSET_KIND),  &
               & 'Viscous diffusion parameter')  )
  CALL pn_check( nfmpi_def_var(ncid_dump, "D_therm", PM_NF_FLOAT,0,0,D_therm_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump,D_therm_varid_dump, 'long_name',  &
               & INT(LEN('Thermal diffusion parameter'),kind=MPI_OFFSET_KIND),  &
               & 'Thermal diffusion parameter')  )
  CALL pn_check( nfmpi_def_var(ncid_dump, "D_comp", PM_NF_FLOAT,0,0,D_comp_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump,D_comp_varid_dump, 'long_name',  &
               & INT(LEN('Compositional diffusion parameter'),kind=MPI_OFFSET_KIND),  &
               & 'Compositional diffusion parameter')  )
  ! PH{
  CALL pn_check( nfmpi_def_var(ncid_dump, "D_mag", PM_NF_FLOAT,0,0,D_mag_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump,D_mag_varid_dump, 'long_name',  &
               & INT(LEN('Magnetic diffusion parameter'),kind=MPI_OFFSET_KIND),  &
               & 'Magnetic diffusion parameter')  )
  ! }PH
  CALL pn_check( nfmpi_def_var(ncid_dump, "S_therm", PM_NF_FLOAT,0,0,S_therm_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump,S_therm_varid_dump, 'long_name',  &
               & INT(LEN('Thermal background stratification parameter'), &
               & kind=MPI_OFFSET_KIND),                                        &
               & 'Thermal background stratification parameter')  )
  CALL pn_check( nfmpi_def_var(ncid_dump, "S_comp", PM_NF_FLOAT,0,0,S_comp_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump,S_comp_varid_dump, 'long_name',  &
               & INT(LEN('Compositional background stratification parameter'), &
               & kind=MPI_OFFSET_KIND),                                        &
               & 'Compositional background stratification parameter')  )
  CALL pn_check( nfmpi_def_var(ncid_dump, "R", PM_NF_FLOAT,0,0,R_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump,R_varid_dump, 'long_name',  &
               & INT(LEN('Rotational parameter'), &
               & kind=MPI_OFFSET_KIND),                                        &
               & 'Rotational parameter')  )
  CALL pn_check( nfmpi_def_var(ncid_dump, "Theta", PM_NF_FLOAT,0,0,Theta_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump,Theta_varid_dump, 'long_name',  &
               & INT(LEN('Angle of rotation axis and gravity parameter'), &
               & kind=MPI_OFFSET_KIND),                                        &
               & 'Angle of rotation axis and gravity parameter')  )
  CALL pn_check( nfmpi_def_var(ncid_dump, "istep", NF_int,0,0,istep_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, istep_varid_dump, 'long_name',&
               & INT(LEN('number of present time step'),kind=MPI_OFFSET_KIND),  &
               & 'number of present time step')  )
  CALL pn_check( nfmpi_def_var(ncid_dump, "time", PM_NF_FLOAT,0,0,time_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, time_varid_dump, 'long_name',&
               & INT(LEN('current model time'),kind=MPI_OFFSET_KIND),       &
               & 'current model time')  )
  CALL pn_check( nfmpi_def_var(ncid_dump, "dt", PM_NF_FLOAT,0,0,dt_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, dt_varid_dump, 'long_name',&
               & INT(LEN('current time step length'),kind=MPI_OFFSET_KIND),    &
               & 'current time step length')  )
  ! Arrays to be stored  
  dimids = (/ ri_dimid_dump, n_dimid_dump, l_dimid_dump, m_dimid_dump /)

  CALL pn_check( nfmpi_def_var(ncid_dump, "ri", NF_int,1,dimids(1),ri_varid_dump) )  
  CALL pn_check( nfmpi_put_att_text(ncid_dump, ri_varid_dump, 'long_name',  &
            & INT(LEN('real and imaginary parts'),kind=MPI_OFFSET_KIND),'real and imaginary parts') ) 
  CALL pn_check( nfmpi_def_var(ncid_dump, "kx", PM_NF_FLOAT,1,dimids(3),kx_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, kx_varid_dump, 'long_name',  &
            & INT(LEN('x-component of wavevector'),kind=MPI_OFFSET_KIND), &
            & 'x-component of wavevector') )
  CALL pn_check( nfmpi_def_var(ncid_dump, "ky", PM_NF_FLOAT,1,dimids(4),ky_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, ky_varid_dump, 'long_name',  &
            & INT(LEN('y-component of wavevector'),kind=MPI_OFFSET_KIND), &
            & 'y-component of wavevector') )
  CALL pn_check( nfmpi_def_var(ncid_dump, "kz", PM_NF_FLOAT,1,dimids(2),kz_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, kz_varid_dump, 'long_name',  &
            & INT(LEN('z-component of wavevector'),kind=MPI_OFFSET_KIND), &
            & 'z-component of wavevector') )

#ifdef TEMPERATURE_FIELD
  CALL pn_check( nfmpi_def_var(ncid_dump, "Temp", PM_NF_FLOAT,4,dimids,Temp_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, Temp_varid_dump, 'long_name',  &
            & INT(LEN('Temperature field'),kind=MPI_OFFSET_KIND),'Temperature field') )
#endif
#ifdef CHEMICAL_FIELD
  CALL pn_check( nfmpi_def_var(ncid_dump, "Chem", PM_NF_FLOAT,4,dimids,Chem_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, Chem_varid_dump, 'long_name',  &
            & INT(LEN('Concentration field'),kind=MPI_OFFSET_KIND),'Concentration field') ) 
#endif

  CALL pn_check( nfmpi_def_var(ncid_dump, "ux", PM_NF_FLOAT,4,dimids,ux_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, ux_varid_dump, 'long_name',  &
            & INT(LEN('x-component of the velocity field' ),kind=MPI_OFFSET_KIND), &
            &'x-component of the velocity field'))
#ifndef TWO_DIMENSIONAL 
  CALL pn_check( nfmpi_def_var(ncid_dump, "uy", PM_NF_FLOAT,4,dimids,uy_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, uy_varid_dump, 'long_name',  &
            & INT(LEN('y-component of the velocity field' ),kind=MPI_OFFSET_KIND), &
            &'y-component of the velocity field'))
#endif 
  CALL pn_check( nfmpi_def_var(ncid_dump, "uz", PM_NF_FLOAT,4,dimids,uz_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, uz_varid_dump, 'long_name',  &
            & INT(LEN('z-component of the velocity field' ),kind=MPI_OFFSET_KIND), &
            &'z-component of the velocity field')) 

#ifdef MAGNETIC
  CALL pn_check( nfmpi_def_var(ncid_dump, "Bx", PM_NF_FLOAT,4,dimids,Bx_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, Bx_varid_dump, 'long_name',  &
            & INT(LEN('x-component of the magnetic field' ),kind=MPI_OFFSET_KIND), &
            &'x-component of the magnetic field'))
#ifndef TWO_DIMENSIONAL 
  CALL pn_check( nfmpi_def_var(ncid_dump, "By", PM_NF_FLOAT,4,dimids,By_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, By_varid_dump, 'long_name',  &
            & INT(LEN('y-component of the magnetic field' ),kind=MPI_OFFSET_KIND), &
            &'y-component of the magnetic field'))
#endif 
  CALL pn_check( nfmpi_def_var(ncid_dump, "Bz", PM_NF_FLOAT,4,dimids,Bz_varid_dump) )
  CALL pn_check( nfmpi_put_att_text(ncid_dump, Bz_varid_dump, 'long_name',  &
            & INT(LEN('z-component of the magnetic field' ),kind=MPI_OFFSET_KIND), &
            &'z-component of the magnetic field'))
#endif
  
  
  ! write global information 
  ! ************************
  CALL pn_check( nfmpi_put_att_text(ncid_dump, nf_global, 'Data',  &
              &  INT(LEN( 'Data needed to restart the DDC code'),kind=MPI_OFFSET_KIND) ,  &
              &  'Data needed to restart the DDC code') )
  IF (myid.EQ.0) THEN
     datestring = ""
     CALL DATE_AND_TIME(real_clock(1),real_clock(2),real_clock(3),date_time)
     WRITE(datestring,'(a,i4,a,i2,a,i2,a,i2,a,i2,a,i2 )')                            &
          &  "Year:",date_time(1),", month:",date_time(2),", day:",date_time(3),     &
          &  ", Time:",date_time(5),":",date_time(6),":",date_time(7)  
  ENDIF
  CALL MPI_BCAST(datestring,LEN(datestring),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  CALL pn_check( nfmpi_put_att_text(ncid_dump, nf_global, 'Created',  &
                                   & INT(LEN(TRIM(datestring)),kind=MPI_OFFSET_KIND) ,TRIM(datestring)) )
  ! switch to data mode and write data to dump file
  ! ***********************************************
  CALL pn_check( nfmpi_enddef(ncid_dump))
  ! write Parameter (only Process 0 needs to do this...)
  CALL pn_check( nfmpi_begin_indep_data(ncid_dump) )
  IF (myid.EQ.0) THEN
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,Gammax_varid_dump,Gammax))
#ifndef TWO_DIMENSIONAL
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,Gammay_varid_dump,Gammay))
#endif
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,Gammaz_varid_dump,Gammaz))
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,B_therm_varid_dump,B_therm))
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,B_comp_varid_dump,B_comp))
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,C_Lorentz_varid_dump,C_Lorentz)) ! PH
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,D_visc_varid_dump,D_visc))
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,D_therm_varid_dump,D_therm))
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,D_comp_varid_dump,D_comp))
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,D_mag_varid_dump,D_mag)) ! PH
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,S_therm_varid_dump,S_therm))
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,S_comp_varid_dump,S_comp))
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,R_varid_dump,R))
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,Theta_varid_dump,Theta))
     CALL pn_check( nfmpi_put_var_int(ncid_dump,ri_varid_dump,(/ (i,i=0,1) /)))
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,kx_varid_dump,kx))
#ifndef TWO_DIMENSIONAL
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,ky_varid_dump,ky))
#endif
     CALL pn_check( PM_NFMPI_PUT_VAR_FLOAT(ncid_dump,kz_varid_dump,kz))
  ENDIF
  CALL pn_check( nfmpi_end_indep_data(ncid_dump) )

  ! flush all buffers 
  CALL pn_check( nfmpi_sync(ncid_dump) )
END SUBROUTINE pn_open_dump
