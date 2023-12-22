FUNCTION d_by_dy(x)
  USE defprecision_module
  USE mpi_transf_module, ONLY:mysx_spec,myex_spec,mysy_spec,myey_spec
  USE parameter_module, ONLY: Lmax,Mmax,Nmax,ky
  IMPLICIT NONE
  COMPLEX(kind=kr) :: d_by_dy(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec)
  COMPLEX(kind=kr) :: x(0:,mysx_spec:,mysy_spec:)
  COMPLEX(kind=kr), PARAMETER :: iu = (0._kr,1._kr) 
  REAL(kind=kr) :: hky
  INTEGER(kind=ki) :: i,j,k,nqfacx,nqfacy,nqfacz,nqfac

#ifdef TWO_DIMENSIONAL
  d_by_dy = 0._kr
#else 
  ! Not cache efficient so far. May need optimization later.
  DO j=mysy_spec,myey_spec
     nqfacy = ( 1 - (j-2*Mmax*((j-1)/(Mmax)))/Mmax )
     hky = ky(j)
     DO i=mysx_spec,myex_spec
        nqfacx = ( 1 - i/Lmax)
        DO k=0,2*Nmax-1
           nqfacz = ( 1 - (k-2*Nmax*((k-1)/(Nmax)))/Nmax )
           nqfac = nqfacx * nqfacy * nqfacz
           d_by_dy(k,i,j) = nqfac * iu * hky * x(k,i,j)
        ENDDO
     ENDDO
  ENDDO
#endif

END FUNCTION d_by_dy
