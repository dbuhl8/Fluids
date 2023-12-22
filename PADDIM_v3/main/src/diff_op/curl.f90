FUNCTION curl(x)
  USE defprecision_module
  USE mpi_transf_module, ONLY:mysx_spec,myex_spec,mysy_spec,myey_spec
  USE parameter_module, ONLY: Gammax,Gammay,Lmax,Mmax,Nmax,kx,ky,kz
  IMPLICIT NONE
  COMPLEX(kind=kr) :: curl(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,dim_curl)
  COMPLEX(kind=kr) :: x(0:,mysx_spec:,mysy_spec:,vec_x:)
  COMPLEX(kind=kr), PARAMETER :: iu = (0._kr,1._kr) 
  REAL(kind=kr) :: hkx,hky,hkz
  INTEGER(kind=ki) :: i,j,k,nqfacx,nqfacy,nqfacz,nqfac

  ! Not cache efficient so far. May need optimization later.

#ifdef TWO_DIMENSIONAL
  j=0
  DO i=mysx_spec,myex_spec
     nqfacx = ( 1 - i/Lmax)
     hkx = kx(i)
     DO k=0,2*Nmax-1
        nqfacz = ( 1 - (k-2*Nmax*((k-1)/(Nmax)))/Nmax )
        hkz = kz(k)
        nqfac = nqfacx * nqfacz
        curl(k,i,j,curl_y) = nqfac * iu * ( hkz * x(k,i,j,vec_x) - hkx * x(k,i,j,vec_z) )
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
           curl(k,i,j,curl_x) = nqfac * iu * ( hky * x(k,i,j,vec_z) - hkz * x(k,i,j,vec_y) )
           curl(k,i,j,curl_y) = nqfac * iu * ( hkz * x(k,i,j,vec_x) - hkx * x(k,i,j,vec_z) )
           curl(k,i,j,curl_z) = nqfac * iu * ( hkx * x(k,i,j,vec_y) - hky * x(k,i,j,vec_x) )
        ENDDO
     ENDDO
  ENDDO
#endif  

END FUNCTION curl

