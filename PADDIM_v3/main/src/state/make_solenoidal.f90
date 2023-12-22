! This subroutines makes a vector field x (given in spectral space) solenoidal. 
! This is done by assuming a correction of type grad(phi), where phi is scalar. 
! Then, the field can be made divergence free by projecting the divergent part out, ie:
!
!   1) solve laplace(phi) = -div(x)
!   2) correct x according to x -> x + grad(phi)
!
! All Nynquist modes are set to zero. 
!
SUBROUTINE make_solenoidal(x)
  USE defprecision_module
  USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec,mysz_spec
  USE parameter_module, ONLY: Lmax,Mmax,Nmax,kx,ky,kz
  IMPLICIT NONE
  COMPLEX(kind=kr) :: x(0:,mysx_spec:,mysy_spec:,vec_x:)
  COMPLEX(kind=kr), PARAMETER :: iu = (0._kr,1._kr) 
  REAL(kind=kr) :: hkx,hky,hkz,ksquare
  INTEGER(kind=ki) :: i,j,k,nqfacx,nqfacy,nqfacz,nqfac
  COMPLEX(kind=kr) :: divu , phi

! Correct velocity field
#ifdef TWO_DIMENSIONAL
! 2D Case 
  j=mysy_spec !=myey_spec=0
  DO i=mysx_spec,myex_spec
     nqfacx = ( 1 - i/Lmax)
     hkx = kx(i)
     DO k=0,2*Nmax-1
        nqfacz = ( 1 - (k-2*Nmax*((k-1)/(Nmax)))/Nmax )
        hkz = kz(k)
        ksquare = hkx**2 + hkz**2
        ksquare =  MAX(ksquare,EPSILON(1._kr)) ! Avoid floating exception in 1/ksquare later..
        nqfac = nqfacx * nqfacz !=0 for Nynquist modes
        ! compute div(x)
        divu = nqfac * iu * ( hkx * x(k,i,j,vec_x) + hkz * x(k,i,j,vec_z) ) 
        ! solve laplace(phi) = -div(x)
        phi = divu / ksquare
        ! correct x according to x -> x + grad(phi)
        x(k,i,j,vec_x) = x(k,i,j,vec_x) + iu * hkx * phi
        x(k,i,j,vec_z) = x(k,i,j,vec_z) + iu * hkz * phi
     ENDDO
  ENDDO
#else 
! 3D Case 
  DO j=mysy_spec,myey_spec
     nqfacy = ( 1 - (j-2*Mmax*((j-1)/(Mmax)))/Mmax )
     hky = ky(j)
     DO i=mysx_spec,myex_spec
        nqfacx = ( 1 - i/Lmax)
        hkx = kx(i)
        DO k=0,2*Nmax-1
           nqfacz = ( 1 - (k-2*Nmax*((k-1)/(Nmax)))/Nmax )
           hkz = kz(k)
           ksquare = hkx**2 + hky**2 + hkz**2
           ksquare =  MAX(ksquare,EPSILON(1._kr)) ! Avoid floating exception in 1/ksquare later..
           nqfac = nqfacx * nqfacy * nqfacz !=0 for Nynquist modes
           ! compute div(x)
           divu = nqfac * iu * ( hkx * x(k,i,j,vec_x) + hky * x(k,i,j,vec_y) + hkz * x(k,i,j,vec_z) ) 
           ! solve laplace(phi) = -div(x)
           phi = divu / ksquare
           ! correct x according to x -> x + grad(phi)
           x(k,i,j,vec_x) = x(k,i,j,vec_x) + iu * hkx * phi
           x(k,i,j,vec_y) = x(k,i,j,vec_y) + iu * hky * phi
           x(k,i,j,vec_z) = x(k,i,j,vec_z) + iu * hkz * phi
        ENDDO
     ENDDO
  ENDDO
#endif

END SUBROUTINE make_solenoidal
