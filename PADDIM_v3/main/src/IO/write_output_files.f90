SUBROUTINE write_output_files(u,Temp,Chem,B,t,dt,istep)
  USE defprecision_module
  USE state_module, ONLY: velocity,buoyancy, field
  USE pnetCDF_IO_module, ONLY: pn_write_dump,pn_write_step_simdat_file,    &
     & pn_write_xyslice_simdat_file,pn_write_xzslice_simdat_file,pn_write_yzslice_simdat_file
  IMPLICIT NONE
  TYPE(velocity) :: u
  TYPE(field)    :: B
  TYPE(buoyancy) :: Temp,Chem
  REAL(kind=kr) :: t,dt
  INTEGER(kind=ki) :: istep

  IF (MOD(istep,n_comp_diag)==0)             &       !Write out diagnostics
       &   CALL write_diagnostics_file(u,Temp,Chem,B,istep,t,dt) ! PH
  IF (MOD(istep,n_wrt_spec)==0)            &         !Write out horizontal spectra
       &   CALL  write_horizontal_spectra(u,Temp,Chem,B,istep,t) ! PH
  IF (MOD(istep,n_wrt_spec)==0)            &          !Write out vertical spectra
       &   CALL  write_vertical_spectra(u,Temp,Chem,B,istep,t) ! PH
  IF (MOD(istep,n_wrt_spec)==0)            &          !Write out shell spectra
       &   CALL  write_shell_spectra(u,Temp,Chem,B,istep,t) ! PH
  IF (MOD(istep,n_wrt_prof)==0)            &          !Write out vertical profiles
       &   CALL  write_z_profile(u,Temp,Chem,B,istep,t) ! PH
  IF (MOD(istep,n_wrt_prof)==0)            &          !Write out y  profiles
       &   CALL  write_y_profile(u,Temp,Chem,B,istep,t) ! PH
  IF (MOD(istep,n_wrt_jc)==0)                &       !Write compressed file  
       &   CALL write_compressed_file(u,Temp,Chem,istep,t,dt) 
  IF (MOD(istep,n_wrt_dump)==0) then  ! DB          &       !Write restart file (netCDF)
           CALL pn_write_dump(u,Temp,Chem,B,t,dt,istep) ! PH
           !CALL write_forcing_file() ! DB write forcing file
  end if ! DB
  IF (MOD(istep,n_wrt_netCDF)==0)            &       !Write simulation data file (netCDF)
       &   CALL pn_write_step_simdat_file(u,Temp,Chem,B,istep,t,dt) ! PH
#ifndef TWO_DIMENSIONAL 
  IF (MOD(istep,n_wrt_slice)==0) THEN                  !Write slice data file (netCDF)
     CALL pn_write_xyslice_simdat_file(u,Temp,Chem,B,istep,t,dt) 
     CALL pn_write_xzslice_simdat_file(u,Temp,Chem,B,istep,t,dt)
     CALL pn_write_yzslice_simdat_file(u,Temp,Chem,B,istep,t,dt) 
  ENDIF
#endif

END SUBROUTINE write_output_files
