SUBROUTINE compute_phys_space_vars(u,Temp,Chem,B)
  USE defprecision_module
  USE parameter_module, ONLY: Nmax
  USE state_module, ONLY: velocity,buoyancy,field,ltime0,shift_time_pointers ! PH
  USE mpi_transf_module, ONLY : FFT_c2r,mysx_spec,myex_spec,mysy_spec,myey_spec
  USE diff_op_module, ONLY: curl,d_by_dx,d_by_dz
  IMPLICIT NONE
  TYPE(velocity) :: u
  TYPE(field)    :: B ! PH
  TYPE(buoyancy) :: Temp,Chem
  COMPLEX(kind=kr), POINTER :: work(:,:,:,:)

! Transform Temp,Chem,u,B to physical space
#ifdef TEMPERATURE_FIELD
  CALL FFT_c2r(Temp%spec(:,:,:,ltime0),Temp%phys)
#endif

#ifdef CHEMICAL_FIELD
  CALL FFT_c2r(Chem%spec(:,:,:,ltime0),Chem%phys)
#endif

CALL FFT_c2r(u%spec(:,:,:,vec_x,ltime0),u%phys(:,:,:,vec_x))
#ifndef TWO_DIMENSIONAL
  CALL FFT_c2r(u%spec(:,:,:,vec_y,ltime0),u%phys(:,:,:,vec_y))
#endif
  CALL FFT_c2r(u%spec(:,:,:,vec_z,ltime0),u%phys(:,:,:,vec_z))

#ifdef MAGNETIC
  CALL FFT_c2r(B%spec(:,:,:,vec_x,ltime0),B%phys(:,:,:,vec_x))
#ifndef TWO_DIMENSIONAL
  CALL FFT_c2r(B%spec(:,:,:,vec_y,ltime0),B%phys(:,:,:,vec_y))
#endif
  CALL FFT_c2r(B%spec(:,:,:,vec_z,ltime0),B%phys(:,:,:,vec_z))
#endif

  ! compute curl(u) in spectral space and transform to physical space
#ifdef TWO_DIMENSIONAL
  ALLOCATE(work(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,curl_y:curl_y))
  work(:,:,:,curl_y) = d_by_dz(u%spec(:,:,:,vec_x,ltime0)) - d_by_dx(u%spec(:,:,:,vec_z,ltime0))
  CALL FFT_c2r(work(:,:,:,curl_y),u%curl(:,:,:,curl_y))
#else
  ALLOCATE(work(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,curl_x:curl_z))
  work = curl(u%spec(:,:,:,:,ltime0))
  CALL FFT_c2r(work(:,:,:,vec_x),u%curl(:,:,:,curl_x))
  CALL FFT_c2r(work(:,:,:,vec_y),u%curl(:,:,:,curl_y))
  CALL FFT_c2r(work(:,:,:,vec_z),u%curl(:,:,:,curl_z))
#endif

#ifdef MAGNETIC
  ! compute curl(B) in spectral space and transform to physical space
#ifdef TWO_DIMENSIONAL
  work(:,:,:,curl_y) = d_by_dz(B%spec(:,:,:,vec_x,ltime0)) - d_by_dx(B%spec(:,:,:,vec_z,ltime0))
  CALL FFT_c2r(work(:,:,:,curl_y),B%curl(:,:,:,curl_y))
#else
  work = curl(B%spec(:,:,:,:,ltime0))
  CALL FFT_c2r(work(:,:,:,vec_x),B%curl(:,:,:,curl_x))
  CALL FFT_c2r(work(:,:,:,vec_y),B%curl(:,:,:,curl_y))
  CALL FFT_c2r(work(:,:,:,vec_z),B%curl(:,:,:,curl_z))
#endif
#endif

  DEALLOCATE(work)
END SUBROUTINE compute_phys_space_vars
