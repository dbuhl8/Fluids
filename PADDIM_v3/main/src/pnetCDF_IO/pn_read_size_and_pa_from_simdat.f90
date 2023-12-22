! This subroutine reads some header information from the output netCDF file
! and determines if everything is consistent with the user given information
! It is usually not recommanded to use the output netCDF file for restarting
! since it is only written in single precision. 
SUBROUTINE pn_read_size_and_pa_from_simdat(error_code)
  USE defprecision_module
  USE parameter_module,   ONLY : Nx,Ny,Nz,Gammax,Gammay,Gammaz,    &
  &                              B_therm,B_comp,C_Lorentz,D_visc,D_therm,D_comp,D_mag,S_therm,S_comp, R, Theta !PH DB
  USE message_passing_module, ONLY: start_mpi,stop_mpi,myid, nprocs1, nprocs2
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  INTEGER (kind=ki) :: Nx_simdat,Ny_simdat,Nz_simdat, nprocs1_simdat, nprocs2_simdat, error_code
  REAL(kind=kr)     :: Gammax_simdat,Gammay_simdat,Gammaz_simdat
  REAL(kind=kr)     :: B_therm_simdat,B_comp_simdat,C_Lorentz_simdat,D_visc_simdat,D_therm_simdat,D_comp_simdat,D_mag_simdat ! PH
  REAL(kind=kr)     :: S_therm_simdat,S_comp_simdat, R_simdat, Theta_simdat
  INTEGER (kind=MPI_OFFSET_KIND) :: idummy
#include "pnetcdf.inc"
!
  ! Open the restart file - read only
  CALL pn_check( nfmpi_open(MPI_COMM_WORLD,TRIM(netCDF_in_simdat_file_name), &
       &  NF_NOWRITE,MPI_INFO_NULL,ncid_simdat                  ) )
  ! get the dimensions
  CALL pn_check (nfmpi_inq_dimid(ncid_simdat,'X',x_dimid_simdat)        )
  CALL pn_check (nfmpi_inq_dimlen(ncid_simdat,x_dimid_simdat,idummy)    )
  Nx_simdat=idummy
  CALL pn_check (nfmpi_inq_dimid(ncid_simdat,'Y',y_dimid_simdat)        )
  CALL pn_check (nfmpi_inq_dimlen(ncid_simdat,y_dimid_simdat,idummy)    )
  Ny_simdat=idummy
  CALL pn_check (nfmpi_inq_dimid(ncid_simdat,'Z',z_dimid_simdat)        )
  CALL pn_check (nfmpi_inq_dimlen(ncid_simdat,z_dimid_simdat,idummy)    )
  Nz_simdat=idummy
  CALL pn_check (nfmpi_inq_dimid(ncid_simdat,'TIME',time_dimid_simdat)     )
  !
  IF (Nx .NE.Nx_simdat  .OR. Ny.NE.Ny_simdat .OR. Nz.NE.Nz_simdat)  THEN
     WRITE(*,'(a,i4,a)') "Simdat file for process ",myid, &
          &                   "is not valid for this problem size."
     WRITE(*,'(a)') "Aborting..."
     CALL stop_mpi
     STOP
  ENDIF

  ! get the varids
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'Gammax',Gammax_varid_simdat)       )
#ifndef TWO_DIMENSIONAL
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'Gammay',Gammay_varid_simdat)       )
#endif
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'Gammaz',Gammaz_varid_simdat)       )
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'B_therm',B_therm_varid_simdat)     )
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'B_comp',B_comp_varid_simdat)       )
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'C_Lorentz',C_Lorentz_varid_simdat) ) ! PH
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'D_visc',D_visc_varid_simdat)       )
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'D_therm',D_therm_varid_simdat)     )
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'D_comp',D_comp_varid_simdat)       )
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'D_mag',D_mag_varid_simdat)         ) ! PH
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'S_therm',S_therm_varid_simdat)     )
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'S_comp',S_comp_varid_simdat)       )
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'R',R_varid_simdat)       )
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'Theta',Theta_varid_simdat)       )
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'t',time_varid_simdat)              )
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'timestep',timestep_varid_simdat)   )
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'dt',dt_varid_simdat)               )
#ifdef TEMPERATURE_FIELD
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'Temp',Temp_varid_simdat)  )
#endif
#ifdef CHEMICAL_FIELD
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'Chem',Chem_varid_simdat)  )
#endif
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'ux',ux_varid_simdat)      )
#ifndef TWO_DIMENSIONAL
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'uy',uy_varid_simdat)      )
#endif
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'uz',uz_varid_simdat)      )
#ifdef MAGNETIC
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'Bx',Bx_varid_simdat)      )
#ifndef TWO_DIMENSIONAL
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'By',By_varid_simdat)      )
#endif
  CALL pn_check (nfmpi_inq_varid(ncid_simdat,'Bz',Bz_varid_simdat)      )
#endif
  
  ! read Parameter values and check if they have changed
  ! write warning to standart output if they have 
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,Gammax_varid_simdat,Gammax_simdat) )
#ifndef TWO_DIMENSIONAL
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,Gammay_varid_simdat,Gammay_simdat) )
#endif
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,Gammaz_varid_simdat,Gammaz_simdat) )

  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,B_therm_varid_simdat,B_therm_simdat) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,B_comp_varid_simdat,B_comp_simdat)   )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,C_Lorentz_varid_simdat,C_Lorentz_simdat)   ) ! PH
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,D_visc_varid_simdat,D_visc_simdat)   )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,D_therm_varid_simdat,D_therm_simdat) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,D_comp_varid_simdat,D_comp_simdat)   )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,D_mag_varid_simdat,D_mag_simdat)     ) ! PH
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,S_therm_varid_simdat,S_therm_simdat) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,S_comp_varid_simdat,S_comp_simdat)   )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,R_varid_simdat,R_simdat)   )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_simdat,Theta_varid_simdat,Theta_simdat)   )
#ifdef STOCHASTIC_FORCING
  CALL pn_check( nfmpi_get_var_int(ncid_simdat, nprocs1_varid_simdat, nprocs1_simdat))
  CALL pn_check( nfmpi_get_var_int(ncid_simdat, nprocs2_varid_simdat, nprocs2_simdat))
#endif 

#ifdef STOCHASTIC_FORCING
  IF (nprocs1_simdat.NE.nprocs1) THEN  !DB
    error_code = 5 !DB
    WRITE(*, '(a, i4, a)') "Warning: Nprocs1 on id ", myid, " is different from value in restart file" !DB
  END IF !DB
  IF (nprocs2_simdat.NE.nprocs2) THEN !DB
    error_code = 6 !DB
    WRITE(*, '(a, i4, a)') "Warning: Nprocs1 on id ", myid, " is different from value in restart file" !DB
  END IF !DB
#endif

  IF (Gammax_simdat.NE.Gammax) WRITE(*,'(a,i4,a,E15.7)') "Warning: Gammax on id ",myid, &
       & " is different from value in restart file. Rel. Error:",                       &
       & ABS((Gammax_simdat-Gammax)/Gammax)
#ifndef TWO_DIMENSIONAL
  IF (Gammay_simdat.NE.Gammay) WRITE(*,'(a,i4,a,E15.7)') "Warning: Gammay on id ",myid, &
       & " is different from value in restart file. Rel. Error:",                       &
       & ABS((Gammay_simdat-Gammay)/Gammay)
#endif
  IF (Gammaz_simdat.NE.Gammaz) WRITE(*,'(a,i4,a,E15.7)') "Warning: Gammaz on id ",myid, &
       & " is different from value in restart file. Rel. Error:",                       &
       & ABS((Gammaz_simdat-Gammaz)/Gammaz)
  IF (B_therm_simdat.NE.B_therm) WRITE(*,'(a,i4,a,E15.7)') "Warning: B_therm on id ",myid, &
       & " is different from value in restart file. Rel. Error:",           &
       & ABS((B_therm_simdat-B_therm)/B_therm)
  IF (B_comp_simdat.NE.B_comp) WRITE(*,'(a,i4,a,E15.7)') "Warning: B_comp on id ",myid, &
       & " is different from value in restart file. Rel. Error:",           &
       & ABS((B_comp_simdat-B_comp)/B_comp)
  IF (C_Lorentz_simdat.NE.C_Lorentz) WRITE(*,'(a,i4,a,E15.7)') "Warning: C_Lorentz on id ",myid, &
       & " is different from value in restart file. Rel. Error:",        &
       & ABS((C_Lorentz_simdat-C_Lorentz)/C_Lorentz) ! PH
  IF (D_visc_simdat.NE.D_visc) WRITE(*,'(a,i4,a,E15.7)') "Warning: D_visc on id ",myid, &
       & " is different from value in restart file. Rel. Error:",        &
       & ABS((D_visc_simdat-D_visc)/D_visc)
  IF (D_therm_simdat.NE.D_therm) WRITE(*,'(a,i4,a,E15.7)') "Warning: D_therm on id ",myid, &
       & " is different from value in restart file. Rel. Error:",        &
       & ABS((D_therm_simdat-D_therm)/D_therm)
  IF (D_comp_simdat.NE.D_comp) WRITE(*,'(a,i4,a,E15.7)') "Warning: D_comp on id ",myid, &
       & " is different from value in restart file. Rel. Error:",        &
       & ABS((D_comp_simdat-D_comp)/D_comp)
  IF (D_mag_simdat.NE.D_mag) WRITE(*,'(a,i4,a,E15.7)') "Warning: D_mag on id ",myid, &
       & " is different from value in restart file. Rel. Error:",        &
       & ABS((D_mag_simdat-D_mag)/D_mag) ! PH
  IF (S_therm_simdat.NE.S_therm) WRITE(*,'(a,i4,a,E15.7)') "Warning: S_therm on id ",myid, &
       & " is different from value in restart file. Rel. Error:",           &
       & ABS((S_therm_simdat-S_therm)/S_therm)
  IF (S_comp_simdat.NE.S_comp) WRITE(*,'(a,i4,a,E15.7)') "Warning: S_comp on id ",myid, &
       & " is different from value in restart file. Rel. Error:",           &
       & ABS((S_comp_simdat-S_comp)/S_comp)
  IF (R_simdat.NE.R) WRITE(*,'(a,i4,a,E15.7)') "Warning: R on id ",myid, &
       & " is different from value in restart file. Rel. Error:",           &
       & ABS((R_simdat-R)/R)
  IF (Theta_simdat.NE.Theta) WRITE(*,'(a,i4,a,E15.7)') "Warning: Theta on id ",myid, &
       & " is different from value in restart file. Rel. Error:",           &
       & ABS((Theta_simdat-Theta)/Theta)

END SUBROUTINE pn_read_size_and_pa_from_simdat
