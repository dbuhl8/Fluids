!+--------------------------------------------------------------------+
!| The following function computes the L2-norm of the velocity,       |
!| vorticity,temperature and chemical field.                          |
!+--------------------------------------------------------------------+
SUBROUTINE compute_uTCB_rms(urms,VORTrms,TEMPrms,CHEMrms,Brms,uxrms,uyrms,uzrms, &
  &                        VORTXrms,VORTYrms,VORTZrms,Bxrms,Byrms,Bzrms,& 
  &                        u,Temp,Chem,B) ! PH
  USE defprecision_module
  USE state_module, ONLY: velocity, buoyancy, field, ltime0 ! PH
  USE message_passing_module, ONLY : myid
  IMPLICIT NONE
  TYPE(velocity)   :: u
  TYPE(field)      :: B
  TYPE(buoyancy)   :: Temp,Chem
  REAL(kind=kr)    :: urms,VORTrms,TEMPrms,CHEMrms, Brms ! PH
  REAL(kind=kr)    :: uxrms,uyrms,uzrms,VORTXrms,VORTYrms,VORTZrms, Bxrms, Byrms, Bzrms ! PH

  ! compute L2-Norm of u
  urms    = rms(u%spec(:,:,:,:,ltime0)) !urms = rms(u%phys) more expensive 
  IF (myid.EQ.0) PRINT*,"kinetic Energy: = ",urms**2

#ifdef MAGNETIC
  ! compute L2-Norm of B
  Brms    = rms(B%spec(:,:,:,:,ltime0))
  IF (myid.EQ.0) PRINT*,"magnetic Energy: = ",Brms**2
#else
  Brms  = 0._kr
#endif

  ! compute L2-Norm of Temperature 
#ifdef TEMPERATURE_FIELD
  TEMPrms = rms(Temp%spec(:,:,:,ltime0)) !  TEMPrms = rms(Temp%phys) more epensive
#else
  TEMPrms = 0._kr
#endif

  ! compute L2-Norm of Chemical field
#ifdef CHEMICAL_FIELD
  CHEMrms = rms(Chem%spec(:,:,:,ltime0)) !  CHEMrms = rms(Chem%phys) more expensive
#else
  CHEMrms = 0._kr
#endif

  ! compute L2-Norm of VORTicity
#ifdef TWO_DIMENSIONAL
  VORTrms = rms_scalar(u%curl(:,:,:,curl_y)) ! No spectral representation stored. Use physical one. 
#else
  VORTrms = rms_vector(u%curl) ! No spectral representation stored. Use physical one. 
#endif

  ! compute L2-Norm of the three components of velocity
  uxrms = rms(u%spec(:,:,:,vec_x,ltime0)) 
#ifndef TWO_DIMENSIONAL
  uyrms = rms(u%spec(:,:,:,vec_y,ltime0)) 
#else
  uyrms = 0._kr
#endif
  uzrms = rms(u%spec(:,:,:,vec_z,ltime0)) 

  ! compute L2-Norm of the three components of vorticity
#ifndef TWO_DIMENSIONAL
  VORTXrms = rms(u%curl(:,:,:,curl_x)) 
  VORTZrms = rms(u%curl(:,:,:,curl_z)) 
#else
  VORTZrms = 0._kr
  VORTXrms = 0._kr
#endif
  VORTYrms = rms(u%curl(:,:,:,curl_y)) 
 
#ifdef MAGNETIC
  ! compute L2-Norm of the three components of the B field
  Bxrms = rms(B%spec(:,:,:,vec_x,ltime0)) 
#ifndef TWO_DIMENSIONAL
  Byrms = rms(B%spec(:,:,:,vec_y,ltime0)) 
#else
  Byrms = 0._kr
#endif
  Bzrms = rms(B%spec(:,:,:,vec_z,ltime0)) 
#else
  Bxrms = 0._kr
  Byrms = 0._kr
  Bzrms = 0._kr
#endif



END SUBROUTINE compute_uTCB_rms
