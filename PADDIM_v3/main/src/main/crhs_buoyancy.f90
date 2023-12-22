! Computes the RHS of the buoyancy equation
SUBROUTINE crhs_buoyancy(rhs,buo,u,ltime_uz,S_buo)
  USE defprecision_module
  USE state_module, ONLY: buoyancy,velocity
  USE parameter_module, ONLY: Nx,Nmax,kx,ky,kz
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys, &
      &                        mysx_spec,myex_spec,mysy_spec,myey_spec, &
      &                        FFT_r2c
  IMPLICIT NONE
  TYPE(buoyancy) :: buo   ! buoyancy field
  TYPE(velocity) :: u     ! velocity field
  COMPLEX(kind=kr) :: rhs(0:,mysx_spec:,mysy_spec:) ! RHS 
  INTEGER(kind=ki) :: ltime_uz ! time level of uz to be used  
  REAL(kind=kr) :: S_buo ! background stratification parameter 
  REAL (kind=kr),POINTER     :: work_phys(:,:,:)
  COMPLEX (kind=kr), POINTER :: work_spec(:,:,:)
  REAL (kind=kr) :: hkx,hky,hkz
  INTEGER (kind=ki) :: i,j,k
  COMPLEX(kind=kr), PARAMETER :: iu = (0._kr,1._kr)

! Compute fourier coefficients of div(buo*u)  

  ALLOCATE(work_phys(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys)) 
  ALLOCATE(work_spec(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec))


  work_phys = buo%phys * u%phys(:,:,:,vec_x)
  CALL FFT_r2c(work_phys,work_spec)
  DO j=mysy_spec,myey_spec
     DO i=mysx_spec,myex_spec
        hkx = kx(i)
        DO k=0,2*Nmax-1
           rhs(k,i,j) = iu * hkx * work_spec(k,i,j) 
        ENDDO
     ENDDO
  ENDDO

#ifndef TWO_DIMENSIONAL
  work_phys = buo%phys * u%phys(:,:,:,vec_y)
  CALL FFT_r2c(work_phys,work_spec)
  DO j=mysy_spec,myey_spec
     hky = ky(j)
     DO i=mysx_spec,myex_spec
        DO k=0,2*Nmax-1
           rhs(k,i,j) = rhs(k,i,j) + iu * hky * work_spec(k,i,j) 
        ENDDO
     ENDDO
  ENDDO
#endif

  work_phys = buo%phys * u%phys(:,:,:,vec_z)
  CALL FFT_r2c(work_phys,work_spec)
  DO j=mysy_spec,myey_spec
     DO i=mysx_spec,myex_spec
        DO k=0,2*Nmax-1
           hkz = kz(k)
           rhs(k,i,j) = rhs(k,i,j) + iu * hkz * work_spec(k,i,j) 
        ENDDO
     ENDDO
  ENDDO

  DEALLOCATE(work_phys,work_spec)

  ! Add term arising from background stratification and change sign 

  rhs = -rhs - S_buo * u%spec(:,:,:,vec_z,ltime_uz)

  ! set to zero constant part 
  IF (mysx_spec.EQ.0 .AND. mysy_spec.EQ.0) rhs(0,0,0) = (0._kr,0._kr)

END SUBROUTINE crhs_buoyancy

