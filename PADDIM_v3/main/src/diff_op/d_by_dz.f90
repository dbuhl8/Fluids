FUNCTION d_by_dz(x)
  USE defprecision_module
  USE mpi_transf_module, ONLY:mysx_spec,myex_spec,mysy_spec,myey_spec
  USE parameter_module, ONLY: Lmax,Mmax,Nmax,kz
  IMPLICIT NONE
  COMPLEX(kind=kr) :: d_by_dz(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec)
  COMPLEX(kind=kr) :: x(0:,mysx_spec:,mysy_spec:)
  COMPLEX(kind=kr), PARAMETER :: iu = (0._kr,1._kr) 
  REAL(kind=kr) :: hkz
  INTEGER(kind=ki) :: i,j,k,nqfacx,nqfacy,nqfacz,nqfac

  ! Not cache efficient so far. May need optimization later.
  DO j=mysy_spec,myey_spec
#ifdef TWO_DIMENSIONAL
     nqfacy = 1
#else
     nqfacy = ( 1 - (j-2*Mmax*((j-1)/(Mmax)))/Mmax )
#endif
     DO i=mysx_spec,myex_spec
        nqfacx = ( 1 - i/Lmax)
        DO k=0,2*Nmax-1
           hkz = kz(k)
           nqfacz = ( 1 - (k-2*Nmax*((k-1)/(Nmax)))/Nmax )
           nqfac = nqfacx * nqfacy * nqfacz
           d_by_dz(k,i,j) = nqfac * iu * hkz * x(k,i,j)
        ENDDO
     ENDDO
  ENDDO
  
END FUNCTION d_by_dz
