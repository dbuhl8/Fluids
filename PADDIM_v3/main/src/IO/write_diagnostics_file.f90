!+--------------------------------------------------------------------+
!| The following subroutine computes various diagnostical parameters  |
!| which are then written to an external output file.                 |
!| Caution: It is assumed that this routine is not called at every    |
!|          time step so that it does not need to be optimized.       |
!|          In fact it is likely to be highly inefficient...          |
!|          To lazy to optimize now, any volunteers welcome!          |
!+--------------------------------------------------------------------+
!| Author: Stephan Stellmach                                          |
!+--------------------------------------------------------------------+
SUBROUTINE write_diagnostics_file(u,Temp,Chem,B,istep,t,dt)
  USE defprecision_module
  USE forcing_module, ONLY: forcing, str
  USE state_module, ONLY: buoyancy,velocity,field,ltime0 ! PH
  USE diagnostics_module, ONLY: compute_uTCB_rms,compute_uTCB_minmax,compute_average_flux, &
      &                         dissipation_buo, compute_power_input, compute_weighted_u_z ! PH, DB
  USE message_passing_module, ONLY : myid
  USE testing_module, ONLY: error_Temp_adv_diff,peak_div_u,peak_div_B ! PH
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys, myex_phys, mysx_phys, &
      &                        mysx_spec,myex_spec,mysy_spec,myey_spec, mysz_spec, myez_spec 
  USE parameter_module, ONLY: Nx,Nmax
  IMPLICIT NONE
  TYPE(velocity)   :: u
  TYPE(field)      :: B
  TYPE(buoyancy)   :: Temp,Chem
  INTEGER(kind=ki) :: istep
  REAL(kind=kr)    :: t,dt
  REAL(kind=kr)    :: urms,VORTrms,TEMPrms,CHEMrms,Brms, flux_Temp,flux_Chem,errU,errB ! PH
  REAL(kind=kr)    :: uxrms,uyrms,uzrms,VORTXrms,VORTYrms,VORTZrms, Bxrms, Byrms, Bzrms ! PH
  REAL(kind=kr)    :: Temp_min,Temp_max,Chem_min,Chem_max,           &
     &                u_min(3),u_max(3),VORT_min(3),VORT_max(3),     &
     &                B_min(3), B_max(3), u_max_abs,VORT_max_abs, B_max_abs ! PH
  REAL(kind=kr)    :: u_z_turb, u_z_lam
  REAL(kind=kr)    :: diss_Temp,diss_Chem, diss_Mom, diss_Mag 
  REAL(kind=kr)    :: powerinput !DB

! compute rms-values
  CALL compute_uTCB_rms(urms,VORTrms,TEMPrms,CHEMrms,Brms,uxrms,uyrms,uzrms, &
  &                     VORTXrms,VORTYrms,VORTZrms,Bxrms,Byrms,Bzrms, &
  &                     u,Temp,Chem,B) ! PH

! compute minimum and maximum values 
  CALL compute_uTCB_minmax(Temp_min,Temp_max,Chem_min,Chem_max,      &
       &                   u_min,u_max,VORT_min,VORT_max,B_min,B_max, &
       &                   u_max_abs,VORT_max_abs,B_max_abs,         &
       &                   u,Temp,Chem,B) ! PH

! compute weighted u_z rms (indicator for small scale turbulence)
  CALL compute_weighted_u_z(u, VORTZrms, u_z_turb, u_z_lam)

! compute fluxes 
#ifdef TEMPERATURE_FIELD
  CALL compute_average_flux(flux_Temp,Temp,u)
#else 
  flux_Temp = 0._kr
#endif
#ifdef STOCHASTIC_FORCING 
  if (t .ne. 0._kr) then
    CALL compute_power_input(powerinput,u) !DB
  else 
    powerinput = 0._kr ! at t = 0 the force_real array doesn't have values in it yet, so this avoids an error at the first timestep
  end if
#else
  powerinput = 0._kr
#endif
#ifdef CHEMICAL_FIELD
  CALL compute_average_flux(flux_Chem,Chem,u)
#else
  flux_Chem = 0._kr
#endif


! compute dissipation rates
#ifdef TEMPERATURE_FIELD
  diss_Temp = dissipation_buo(Temp%spec(:,:,:,ltime0))
#else
  diss_Temp = 0._kr
#endif
#ifdef CHEMICAL_FIELD
  diss_Chem = dissipation_buo(Chem%spec(:,:,:,ltime0))
#else
  diss_Chem = 0._kr
#endif
#ifdef TWO_DIMENSIONAL
diss_Mom = dissipation_buo(u%spec(:,:,:,vec_x,ltime0)) + dissipation_buo(u%spec(:,:,:,vec_z,ltime0))
#else 
diss_Mom = dissipation_buo(u%spec(:,:,:,vec_x,ltime0)) + dissipation_buo(u%spec(:,:,:,vec_y,ltime0)) &
         + dissipation_buo(u%spec(:,:,:,vec_z,ltime0))
#endif
#ifdef MAGNETIC
#ifdef TWO_DIMENSIONAL
diss_Mag = dissipation_buo(B%spec(:,:,:,vec_x,ltime0)) + dissipation_buo(u%spec(:,:,:,vec_z,ltime0))
#else
diss_Mag = dissipation_buo(B%spec(:,:,:,vec_x,ltime0)) + dissipation_buo(B%spec(:,:,:,vec_y,ltime0)) &
         + dissipation_buo(B%spec(:,:,:,vec_z,ltime0))
#endif
#else
diss_Mag = 0._kr
#endif


! check for peak in spectrum of div(u)

  !call error_Temp_adv_diff(err,Temp,t)
  CALL peak_div_u(u,errU)
#ifdef MAGNETIC
  CALL peak_div_B(B,errB) ! PH
#endif
  

  ! PH{
  IF (myid.EQ.0) THEN
     PRINT*,"Peak spectrum div(u) = ",errU
#ifdef MAGNETIC
     PRINT*,"Peak spectrum div(B) = ",errB
#endif
     WRITE(uout(1),'(i7,50E20.7)')                                                            &
          !         1     2 3  4    5       6       7       8       9         10
          &         istep,t,dt,urms,VORTrms,TEMPrms,CHEMrms,Brms, flux_Temp,flux_Chem,        &
          !         11       12       13       14 
          &         Temp_min,Temp_max,Chem_min,Chem_max,                                      &
          !         15       16       17       18       19       20 
          &         u_min(1),u_max(1),u_min(2),u_max(2),u_min(3),u_max(3),                    &
          !         21           22           23           24            25         26
          &         VORT_min(1),VORT_max(1),VORT_min(2),VORT_max(2), VORT_min(3),VORT_max(3), &
          !         27          28        29        30        31        32
          &         B_min(1), B_max(1), B_min(2), B_max(2), B_min(3), B_max(3),               &
          !         33          34             35
          &         u_max_abs, VORT_max_abs, B_max_abs,                                       &
          !         36     37    38      39      40       41       42     43     44
          &         uxrms,uyrms,uzrms,VORTXrms,VORTYrms,VORTZrms, Bxrms, Byrms, Bzrms,        &
          !         45          46        47       48         49         50       51
          &         diss_Temp,diss_Chem,diss_Mom,diss_Mag,powerinput, u_z_turb, u_z_lam
  ENDIF
  ! }PH

  IF (myid.EQ.0) THEN                  ! Use for testing routines which compute power spectra
     PRINT*,"urms**2=",urms**2
     PRINT*,"Temp_rms**2=",TEMPrms**2
     PRINT*,"Chem_rms**2_",CHEMrms**2
  ENDIF

END SUBROUTINE write_diagnostics_file
