! Compute the RHS of the magnetic equation
SUBROUTINE crhs_magnetic(rhs,u,B)
  USE defprecision_module
  USE state_module, ONLY: buoyancy,velocity, field
  USE parameter_module, ONLY: Nx,Nmax,kx,ky,kz,B_therm,B_comp
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys, &
      &                        mysx_spec,myex_spec,mysy_spec,myey_spec, &
      &                        FFT_r2c
  IMPLICIT NONE
  TYPE(buoyancy) :: Temp,Chem
  TYPE(velocity) :: u
  TYPE(field) :: B
  COMPLEX(kind=kr) ::  rhs(0:,mysx_spec:,mysy_spec:,1:) ! RHS 
  REAL (kind=kr),POINTER     :: work_phys_x(:,:,:), work_phys_y(:,:,:), work_phys_z(:,:,:)
  COMPLEX (kind=kr),POINTER  :: work_spec_x(:,:,:), work_spec_y(:,:,:), work_spec_z(:,:,:)
  COMPLEX(kind=kr), PARAMETER :: iu = (0._kr,1._kr) 
  REAL (kind=kr) :: hkx,hky,hkz
  REAL (kind=kr) :: ksquare
  COMPLEX (kind=kr) :: kN,fac
  INTEGER (kind=ki) :: i,j,k

  ! compute necessary u x B components and transform to spectral space
  ! (for 2D case, "y" component of u x B is the only non-zero one)
  ALLOCATE(work_phys_y(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys))
  ALLOCATE(work_spec_y(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
  work_phys_y = u%phys(:,:,:,vec_z) * B%phys(:,:,:,vec_x) &
            & - u%phys(:,:,:,vec_x) * B%phys(:,:,:,vec_z) 
  CALL FFT_r2c(work_phys_y,work_spec_y)

#ifndef TWO_DIMENSIONAL
  ALLOCATE(work_phys_x(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys))
  ALLOCATE(work_phys_z(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys))
  ALLOCATE(work_spec_x(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
  ALLOCATE(work_spec_z(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
  work_phys_x = u%phys(:,:,:,vec_y) * B%phys(:,:,:,vec_z) &
            & - u%phys(:,:,:,vec_z) * B%phys(:,:,:,vec_y) 
  work_phys_z = u%phys(:,:,:,vec_x) * B%phys(:,:,:,vec_y) &
            & - u%phys(:,:,:,vec_y) * B%phys(:,:,:,vec_x) 
  CALL FFT_r2c(work_phys_x,work_spec_x)
  CALL FFT_r2c(work_phys_z,work_spec_z)
#endif
  
  ! curl u x B in spectral space and store it in rhs array
#ifndef TWO_DIMENSIONAL
  DO j=mysy_spec,myey_spec
     DO i=mysx_spec,myex_spec
        DO k=0,2*Nmax-1
           rhs(k,i,j,vec_x) =  iu*ky(j)*work_spec_z(k,i,j) - iu*kz(k)*work_spec_y(k,i,j)
           rhs(k,i,j,vec_y) =  iu*kz(k)*work_spec_x(k,i,j) - iu*kx(i)*work_spec_z(k,i,j)
           rhs(k,i,j,vec_z) =  iu*kx(i)*work_spec_y(k,i,j) - iu*ky(j)*work_spec_x(k,i,j)
        ENDDO
     ENDDO
  ENDDO
  DEALLOCATE(work_phys_x, work_phys_y, work_phys_z, work_spec_x, work_spec_y, work_spec_z)
#else
  j=mysy_spec !=myey_spec=0
  DO i=mysx_spec,myex_spec
     DO k=0,2*Nmax-1
        rhs(k,i,j,vec_x) = - iu*kz(k)*work_spec_y(k,i,j)
        rhs(k,i,j,vec_z) =   iu*kx(i)*work_spec_y(k,i,j)
     ENDDO
  ENDDO
  DEALLOCATE(work_phys_y, work_spec_y)
#endif

END SUBROUTINE crhs_magnetic
