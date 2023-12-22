! Given the velocity field u, this subroutine computes an estimate for the CFL-Limit
! for the third order Adams-Bashforth/Backward-Differencing (AB/BDF3) scheme.
! The estimation is based on results for the advection-diffusion equation given in 
! the book of Peyret(2002). 
FUNCTION CFL_AB_BDF3(u)
  USE defprecision_module
  USE state_module, ONLY: velocity
  USE parameter_module, ONLY: Nx,Nmax,D_visc,D_therm,D_comp,D_mag,kx,ky,kz ! PH
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys, &
      &                        mysx_spec,myex_spec,mysy_spec,myey_spec
  USE MPI
  IMPLICIT NONE
  TYPE(velocity) :: u
  REAL(kind=kr) :: CFL_AB_BDF3
  REAL(kind=kr) :: lumax(dim_vec),gumax(dim_vec),ku,ksquared,la(4),ga(4)
  REAL(kind=kr) :: hkx,hky,hkz
  REAL(kind=kr), PARAMETER :: a=5.61_kr, b=0.69
  INTEGER(kind=ki) :: l,i,j,k,ierr
  ! Find the maximm of ux,uy and uz
  lumax = 0._kr
  DO l=vec_x,vec_z
     DO k=mysz_phys,myez_phys
        DO j=mysy_phys,myey_phys
           DO i=0,Nx-1
              lumax(l) = MAX(ABS(u%phys(i,j,k,l)),lumax(l))
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  CALL MPI_ALLREDUCE(lumax,gumax,dim_vec,PM_MPI_FLOAT_TYPE,MPI_MAX,MPI_COMM_WORLD,ierr)
  
  ! The stability region boundary is nearly rectilinear and may be parameterized according
  ! to Peyret (2002). 
  la = 0._kr
  DO j=mysy_spec,myey_spec
     hky = ky(j)
     DO i=mysx_spec,myex_spec
        hkx = kx(i)
        DO k=0,2*Nmax-1
           hkz = kz(k)
#ifdef TWO_DIMENSIONAL
           ku = ABS(hkx * gumax(vec_x) + hkz*gumax(vec_z))
           ksquared = hkx**2 + hkz**2
#else
           ku = ABS(hkx * gumax(vec_x) + hky*gumax(vec_y) + hkz*gumax(vec_z))
           ksquared = hkx**2 + hky**2 + hkz**2
#endif
           la(1) = MAX(la(1), a * ku - D_visc *ksquared      ) ! velocity
#ifdef TEMPERATURE_FIELD
           la(2) = MAX(la(2), a * ku - D_therm*ksquared      ) ! temperature
#else 
           la(2) = 0._kr
#endif
#ifdef CHEMICAL_FIELD
           la(3) = MAX(la(3), a * ku - D_comp*ksquared       ) ! concentration
#else
           la(3) = 0._kr
#endif
#ifdef MAGNETIC
           la(4) = MAX(la(4), a * ku - D_mag*ksquared       ) ! magnetic (PH)
#else
           la(4) = 0._kr
#endif
        ENDDO
     ENDDO
  ENDDO
  CALL MPI_ALLREDUCE(la,ga,4,PM_MPI_FLOAT_TYPE,MPI_MAX,MPI_COMM_WORLD,ierr)
  CFL_AB_BDF3 = a * b / MAX(MAXVAL(ga),EPSILON(1._kr))
 

  ! The scheme may be unconditionally stable. In this case, return at huge value
  IF (CFL_AB_BDF3.LE.0._kr) CFL_AB_BDF3=1.E20_kr

END FUNCTION CFL_AB_BDF3
