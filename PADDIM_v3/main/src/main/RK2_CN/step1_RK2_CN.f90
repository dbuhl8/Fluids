! Ths subroutine computes 
! phi1 = (1 - dt * k**2 * dcoeff) phi + dt * rhs
SUBROUTINE step1_RK2_CN(phi1,phi,rhs,dt,dcoeff)
  USE defprecision_module
  USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec, &
      &                        mysy_phys,myey_phys,mysz_phys,myez_phys
  USE parameter_module, ONLY : Lmax,Mmax,Nmax,Nx,kx,ky,kz
  IMPLICIT NONE
  COMPLEX(kind=kr) :: phi1(0:,mysx_spec:,mysy_spec:)
  COMPLEX(kind=kr) ::  phi(0:,mysx_spec:,mysy_spec:)
  COMPLEX(kind=kr) ::  rhs(0:,mysx_spec:,mysy_spec:)
  REAL(kind=kr) :: dt,dcoeff
  REAL(kind=kr) :: fac,kxsquared,kysquared,kzsquared
  INTEGER(kind=ki) :: i,j,k,nqfac,nqfacx,nqfacy,nqfacz

#ifdef TWO_DIMENSIONAL
  j=mysy_spec !=myey_spec=0
  DO i=mysx_spec,myex_spec
     nqfacx = ( 1 - i/Lmax)
     kxsquared = kx(i)**2
     DO k=0,2*Nmax-1
        nqfacz = ( 1 - (k-2*Nmax*((k-1)/(Nmax)))/Nmax )
        kzsquared = kz(k)**2
        nqfac = nqfacx * nqfacz
        fac = ( 1._kr - dt * dcoeff * ( kxsquared + kzsquared ) )
        phi1(k,i,j) = fac * phi(k,i,j) + dt * rhs(k,i,j)
     ENDDO
  ENDDO
#else
  DO j=mysy_spec,myey_spec
     nqfacy = ( 1 - (j-2*Mmax*((j-1)/(Mmax)))/Mmax )
     kysquared = ky(j)**2
     DO i=mysx_spec,myex_spec
        nqfacx = ( 1 - i/Lmax)
        kxsquared = kx(i)**2
        DO k=0,2*Nmax-1
           nqfacz = ( 1 - (k-2*Nmax*((k-1)/(Nmax)))/Nmax )
           kzsquared = kz(k)**2
           nqfac = nqfacx * nqfacy * nqfacz
           !This is to check that the nqfac stuff works...
           !if (nqfac.eq.0. .and. (i.ne.Lmax .and. j.ne.Mmax .and. k.ne.Nmax)) stop
           !if (nqfac.ne.0. .and. (i.eq.Lmax .or. j.eq.Mmax .or. k.eq.Nmax .or. nqfac.ne.1)) stop
           fac = ( 1._kr - dt * dcoeff * ( kxsquared + kysquared + kzsquared ) )
           phi1(k,i,j) = fac * phi(k,i,j) + dt * rhs(k,i,j)
        ENDDO
     ENDDO
  ENDDO
#endif

END SUBROUTINE step1_RK2_CN
