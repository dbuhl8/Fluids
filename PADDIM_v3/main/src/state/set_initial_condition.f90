SUBROUTINE set_initial_condition(u,Temp,Chem,B)
  USE defprecision_module
  USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec,mysz_spec,FFT_r2c,FFT_c2r
  USE parameter_module, ONLY: Nmax 
  USE diff_op_module, ONLY: curl
  IMPLICIT NONE
  TYPE(velocity) :: u
  TYPE(field)    :: B
  TYPE(buoyancy) :: Temp,Chem
  COMPLEX(kind=kr), POINTER :: work(:,:,:,:)
  
! transform to spectral space 

#ifdef TEMPERATURE_FIELD
  Temp%spec=(0._kr,0._kr)
  CALL FFT_r2c(Temp%phys,Temp%spec(:,:,:,ltime0))
#endif

#ifdef CHEMICAL_FIELD
  Chem%spec=(0._kr,0._kr)
  CALL FFT_r2c(Chem%phys,Chem%spec(:,:,:,ltime0))
#endif

  u%spec=(0._kr,0._kr)
#ifdef TWO_DIMENSIONAL
  CALL FFT_r2c(u%phys(:,:,:,vec_x),u%spec(:,:,:,vec_x,ltime0))
  CALL FFT_r2c(u%phys(:,:,:,vec_z),u%spec(:,:,:,vec_z,ltime0))
#else 
  CALL FFT_r2c(u%phys(:,:,:,vec_x),u%spec(:,:,:,vec_x,ltime0))
  CALL FFT_r2c(u%phys(:,:,:,vec_y),u%spec(:,:,:,vec_y,ltime0))
  CALL FFT_r2c(u%phys(:,:,:,vec_z),u%spec(:,:,:,vec_z,ltime0))
#endif
  
#ifdef MAGNETIC
  B%spec=(0._kr,0._kr)
#ifdef TWO_DIMENSIONAL
  CALL FFT_r2c(B%phys(:,:,:,vec_x),B%spec(:,:,:,vec_x,ltime0)) 
  CALL FFT_r2c(B%phys(:,:,:,vec_z),B%spec(:,:,:,vec_z,ltime0)) 
#else
  CALL FFT_r2c(B%phys(:,:,:,vec_x),B%spec(:,:,:,vec_x,ltime0)) 
  CALL FFT_r2c(B%phys(:,:,:,vec_y),B%spec(:,:,:,vec_y,ltime0))
  CALL FFT_r2c(B%phys(:,:,:,vec_z),B%spec(:,:,:,vec_z,ltime0))
#endif
#endif  

! Make initial velocity solenoidal and recompute physical space values
  CALL make_solenoidal(u%spec(:,:,:,:,ltime0))
  CALL FFT_c2r(u%spec(:,:,:,vec_x,ltime0),u%phys(:,:,:,vec_x))
  CALL FFT_c2r(u%spec(:,:,:,vec_z,ltime0),u%phys(:,:,:,vec_z))
#ifndef TWO_DIMENSIONAL
  CALL FFT_c2r(u%spec(:,:,:,vec_y,ltime0),u%phys(:,:,:,vec_y))
#endif

#ifdef MAGNETIC
! Make initial field solenoidal and recompute physical space values
  CALL make_solenoidal(B%spec(:,:,:,:,ltime0))
  CALL FFT_c2r(B%spec(:,:,:,vec_x,ltime0),B%phys(:,:,:,vec_x)) 
  CALL FFT_c2r(B%spec(:,:,:,vec_z,ltime0),B%phys(:,:,:,vec_z)) 
#ifndef TWO_DIMENSIONAL
  CALL FFT_c2r(B%spec(:,:,:,vec_y,ltime0),B%phys(:,:,:,vec_y)) 
#endif
#endif

! compute curl(u) in spectral space and transform to physical space
  ALLOCATE(work(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,dim_curl))
  work = curl(u%spec(:,:,:,:,ltime0))
  CALL FFT_c2r(work(:,:,:,curl_y),u%curl(:,:,:,curl_y))
#ifndef TWO_DIMENSIONAL
  CALL FFT_c2r(work(:,:,:,curl_x),u%curl(:,:,:,curl_x))
  CALL FFT_c2r(work(:,:,:,curl_z),u%curl(:,:,:,curl_z))
#endif
  DEALLOCATE(work)

#ifdef MAGNETIC
! compute curl(B) in spectral space and transform to physical space
  ALLOCATE(work(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,dim_curl))
  work = curl(B%spec(:,:,:,:,ltime0))
  CALL FFT_c2r(work(:,:,:,curl_y),B%curl(:,:,:,curl_y))
#ifndef TWO_DIMENSIONAL
  CALL FFT_c2r(work(:,:,:,curl_x),B%curl(:,:,:,curl_x))
  CALL FFT_c2r(work(:,:,:,curl_z),B%curl(:,:,:,curl_z))
#endif
  DEALLOCATE(work)
#endif
 
! set to zero mean modes (except for MF? ) 
  IF (mysx_spec.EQ.0 .AND. mysy_spec.EQ.0 .AND. mysz_spec.EQ.0) THEN 
     u%spec(0,0,0,:,:) = (0._kr,0._kr)
#ifdef TEMPERATURE_FIELD
     Temp%spec(0,0,0,:) = (0._kr,0._kr)
#endif
#ifdef CHEMICAL_FIELD
     Chem%spec(0,0,0,:) = (0._kr,0._kr)
#endif
  ENDIF

#ifdef AB_BDF3
! set right hand sides to zero
  u%rhs = (0._kr,0._kr)
#ifdef MAGNETIC
  B%rhs = (0._kr,0._kr) 
#endif
#ifdef TEMPERATURE_FIELD
  Temp%rhs = (0._kr,0._kr)
#endif
#ifdef CHEMICAL_FIELD
  Chem%rhs = (0._kr,0._kr)
#endif
#endif

END SUBROUTINE set_initial_condition
