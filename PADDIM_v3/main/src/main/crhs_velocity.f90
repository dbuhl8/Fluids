! Compute the RHS of the velocity equation
SUBROUTINE crhs_velocity(t,rhs,u,Temp,Chem,B)
  USE defprecision_module
  USE defs_2D_3D_module !DB 
  USE forcing_module ! DB
  USE message_passing_module, ONLY: myid ! DB
  USE state_module, ONLY: buoyancy,velocity,field ! PH
  USE parameter_module, ONLY: Nx,Nmax,kx,ky,kz,B_therm,B_comp,C_Lorentz, R, Theta, Gammaz, Gammay, Ny, Nz, pi
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys, myex_phys, mysx_phys, &
      &                        mysx_spec,myex_spec,mysy_spec,myey_spec, mysz_spec, myez_spec, &
      &                        FFT_r2c
  IMPLICIT NONE
  TYPE(buoyancy) :: Temp,Chem
  TYPE(velocity) :: u
  TYPE(field) :: B
  COMPLEX(kind=kr) ::  rhs(0:,mysx_spec:,mysy_spec:,1:)!, force_spec(0:,mysx_spec:,mysy_spec:,1:)
  COMPLEX(kind=kr),POINTER :: force_spec(:,:,:,:)
  REAL (kind=kr),POINTER     :: work_phys(:,:,:) 
  REAL (kind=kr) :: hkx,hky,hkz, dz,zc
  REAL (kind=kr) :: ksquare
  COMPLEX (kind=kr) :: kN,fac
  INTEGER (kind=ki) :: i,j,k
  REAL(kind=kr) :: yc,dy
  REAL(kind=kr) :: t

print *, "cpu "//trim(str(myid))//": about to go through crhs_velocity"

! compute Fourier coeff. of -curl(u) \times u and curl(B) \times B, add buoyancy force if present
  ALLOCATE(work_phys(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys))
  !ALLOCATE(force_spec(0:,mysx_spec:,mysy_spec:,1:))
  ALLOCATE(force_spec(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,vec_x:vec_z))
print *, "cpu "//trim(str(myid))//": about to call forcing.f90"
  CALL forcing(t, force_spec) !DB if stochastic forcing, returns gp forcing, if not returns sin forcing 

print *, "cpu "//trim(str(myid))//":  mysx:", mysx_spec, ", myex:", myex_spec, ", mysy:", mysy_spec, &
        & " myey_spec:", myey_spec, ", mysz:", mysz_spec, " myez_spec:", myez_spec, ", 2*nmax-1:, ", 2*Nmax-1



!print *, "cpu "//trim(str(myid))//": called forcing.f90"

#ifdef TWO_DIMENSIONAL
  ! x-component
  work_phys = - u%curl(:,:,:,curl_y) * u%phys(:,:,:,vec_z) 
#ifdef MAGNETIC
  work_phys = work_phys + C_Lorentz*(B%curl(:,:,:,curl_y) * B%phys(:,:,:,vec_z) )
#endif
  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_x))
  ! z-component 
  work_phys =   u%curl(:,:,:,curl_y) * u%phys(:,:,:,vec_x) 
#ifdef MAGNETIC
  work_phys = work_phys - C_Lorentz*(B%curl(:,:,:,curl_y) * B%phys(:,:,:,vec_x) )
#endif

#else
  ! x-component
  work_phys = - u%curl(:,:,:,curl_y) * u%phys(:,:,:,vec_z) &
            & + u%curl(:,:,:,curl_z) * u%phys(:,:,:,vec_y) &
            & - R * u%phys(:,:,:,vec_z) * sin(Theta) &
            & + R * u%phys(:,:,:,vec_y) * cos(Theta) 
#ifdef MAGNETIC
  work_phys = work_phys + C_Lorentz*(B%curl(:,:,:,curl_y) * B%phys(:,:,:,vec_z) &
            &              - B%curl(:,:,:,curl_z) * B%phys(:,:,:,vec_y) )  
#endif
  
  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_x))

  ! y-component 
  work_phys = - u%curl(:,:,:,curl_z) * u%phys(:,:,:,vec_x) &
            & + u%curl(:,:,:,curl_x) * u%phys(:,:,:,vec_z) &
            & - R * u%phys(:,:,:,vec_x) * cos(Theta)
#ifdef MAGNETIC
  work_phys = work_phys + C_Lorentz*(B%curl(:,:,:,curl_z) * B%phys(:,:,:,vec_x) &
            &              - B%curl(:,:,:,curl_x) * B%phys(:,:,:,vec_z) )
#endif
  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_y))
  ! z-component 
  work_phys = - u%curl(:,:,:,curl_x) * u%phys(:,:,:,vec_y) &
            & + u%curl(:,:,:,curl_y) * u%phys(:,:,:,vec_x) &
            & + R * u%phys(:,:,:,vec_x) * sin(Theta)
#ifdef MAGNETIC
  work_phys = work_phys + C_Lorentz*(B%curl(:,:,:,curl_x) * B%phys(:,:,:,vec_y) &
            &              - B%curl(:,:,:,curl_y) * B%phys(:,:,:,vec_x) )
#endif
#endif
 


#ifdef TEMPERATURE_FIELD
  work_phys = work_phys + B_therm * Temp%phys 
#endif
#ifdef CHEMICAL_FIELD
  work_phys = work_phys - B_comp * Chem%phys
#endif

  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_z))

print *, "cpu "//trim(str(myid))//": about to add force_spec to rhs"
  rhs = rhs + force_spec 
print *, "cpu "//trim(str(myid))//": added force_spec to rhs"
  DEALLOCATE(work_phys)
print *, "cpu "//trim(str(myid))//": added work_phys deallocated"
 DEALLOCATE(force_spec)
print *, "cpu "//trim(str(myid))//": deallocated force_spec"

print *, "cpu "//trim(str(myid))//": got through crhs_velocity"

! add contribution from pressure gradient 
  ! This does not look too cache efficient...
#ifdef TWO_DIMENSIONAL
  j=mysy_spec !=myey_spec=0
  DO i=mysx_spec,myex_spec
     hkx = kx(i)
     DO k=0,2*Nmax-1
        hkz = kz(k)
        kN =   hkx * rhs(k,i,j,vec_x) &
             + hkz * rhs(k,i,j,vec_z)
        ksquare = hkx**2 + hkz**2
        ksquare = MAX(ksquare,EPSILON(1._kr)) ! Avoid floating exception in 1/ksquare later...
        fac = kN / ksquare
        rhs(k,i,j,vec_x) =   rhs(k,i,j,vec_x) - hkx * fac
        rhs(k,i,j,vec_z) =   rhs(k,i,j,vec_z) - hkz * fac
     ENDDO
  ENDDO
#else
  DO j=mysy_spec,myey_spec
     hky = ky(j)
     DO i=mysx_spec,myex_spec
        hkx = kx(i)
        DO k=0,2*Nmax-1
           hkz = kz(k)
           kN =   hkx * rhs(k,i,j,vec_x) &
           &    + hky * rhs(k,i,j,vec_y) &
           &    + hkz * rhs(k,i,j,vec_z)
           ksquare = hkx**2 + hky**2 + hkz**2
           ksquare = MAX(ksquare,EPSILON(1._kr)) ! Avoid floating exception in 1/ksquare later...
           fac = kN / ksquare
           rhs(k,i,j,vec_x) =   rhs(k,i,j,vec_x) - hkx * fac
           rhs(k,i,j,vec_y) =   rhs(k,i,j,vec_y) - hky * fac
           rhs(k,i,j,vec_z) =   rhs(k,i,j,vec_z) - hkz * fac
        ENDDO
     ENDDO
  ENDDO
#endif
  ! Set to zero constant part. 
  IF (mysx_spec.EQ.0 .AND. mysy_spec.EQ.0) rhs(0,0,0,:) = (0._kr,0._kr)

END SUBROUTINE crhs_velocity
