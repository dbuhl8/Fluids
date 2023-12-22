FUNCTION div(x)
  USE defprecision_module
  USE mpi_transf_module, ONLY:mysx_spec,myex_spec,mysy_spec,myey_spec
  USE parameter_module, ONLY: Lmax,Mmax,Nmax,kx,ky,kz
  IMPLICIT NONE
  COMPLEX(kind=kr) :: div(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec)
  COMPLEX(kind=kr) :: x(0:,mysx_spec:,mysy_spec:,1:)
  COMPLEX(kind=kr), PARAMETER :: iu = (0._kr,1._kr) 
  REAL(kind=kr) :: hkx,hky,hkz
  INTEGER(kind=ki) :: i,j,k,nqfacx,nqfacy,nqfacz,nqfac

#ifdef TWO_DIMENSIONAL
  j=0
  DO i=mysx_spec,myex_spec
     nqfacx = ( 1 - i/Lmax)
     hkx = kx(i)
     DO k=0,2*Nmax-1
        nqfacz = ( 1 - (k-2*Nmax*((k-1)/(Nmax)))/Nmax )
        hkz = kz(k)
        nqfac = nqfacx * nqfacz
           div(k,i,j) = nqfac * iu * ( hkx * x(k,i,j,vec_x) + hkz * x(k,i,j,vec_z) ) 
        ENDDO
     ENDDO
#else
  DO j=mysy_spec,myey_spec
     nqfacy = ( 1 - (j-2*Mmax*((j-1)/(Mmax)))/Mmax )
     hky = ky(j)
     DO i=mysx_spec,myex_spec
        nqfacx = ( 1 - i/Lmax)
        hkx = kx(i)
        DO k=0,2*Nmax-1
           nqfacz = ( 1 - (k-2*Nmax*((k-1)/(Nmax)))/Nmax )
           hkz = kz(k)
           nqfac = nqfacx * nqfacy * nqfacz
           div(k,i,j) = nqfac * iu * ( hkx * x(k,i,j,vec_x)  &
           &                         + hky * x(k,i,j,vec_y)  & 
           &                         + hkz * x(k,i,j,vec_z) ) 
        ENDDO
     ENDDO
  ENDDO
#endif

END FUNCTION div
