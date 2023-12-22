! Ths subroutine computes 
! phi = 0.5 * ( (phi + phi1) + dt * rhs ) / (1 + 0.5 * dt * k**2 * dcoeff) 
SUBROUTINE step2_RK2_CN(phi,phin,rhs,dt,dcoeff)
  USE defprecision_module
  USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec, &
      &                        mysy_phys,myey_phys,mysz_phys,myez_phys
  USE parameter_module, ONLY : Lmax,Mmax,Nmax,Nx,kx,ky,kz
  IMPLICIT NONE
  COMPLEX(kind=kr) ::  phi(0:,mysx_spec:,mysy_spec:)
  COMPLEX(kind=kr) :: phin(0:,mysx_spec:,mysy_spec:)
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
        fac = 0.5_kr / ( 1._kr + 0.5 * dt * dcoeff * ( kxsquared + kzsquared ) )
        phi(k,i,j) = nqfac * ( fac * ((phi(k,i,j) + phin(k,i,j)) + dt * rhs(k,i,j)) )
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
           fac = 0.5_kr / ( 1._kr + 0.5 * dt * dcoeff * ( kxsquared + kysquared + kzsquared ) )
           phi(k,i,j) = nqfac * (fac * ((phi(k,i,j) + phin(k,i,j)) + dt * rhs(k,i,j)) )
        ENDDO
     ENDDO
  ENDDO
#endif

END SUBROUTINE step2_RK2_CN
