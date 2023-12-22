SUBROUTINE init_fft_storage_scheme
!+--------------------------------------------------------------------+
!| This subroutine sets up the matrices which contain the FFT storage |
!| scheme.                                                            |
!+--------------------------------------------------------------------+
  IMPLICIT NONE
  INTEGER (kind=ki) :: i,j,k

  ! kx(i) containes the x-wavenumber for modes which are stored at x(:,i,...)
  DO i=0,Lmax
     kx(i)=i*alpha
  ENDDO
  ! ky(j) containes the y-wavenumber for modes which are stored at x(:,:,j,...)
  IF (Mmax.eq.0) then 
     ky = 0._kr 
  ELSE 
     DO j=0,2*Mmax-1
        ky(j)=(j-2*Mmax*((j-1)/(Mmax)))*beta
     ENDDO
  ENDIF
  ! kz(k) containes the z-wavenumber for modes which are stored at (k,:,:,...)
  IF (Nmax.eq.0) then 
     kz = 0._kr 
  ELSE
     DO k=0,2*Nmax-1
        kz(k)=(k-2*Nmax*((k-1)/(Nmax)))*gamma
     ENDDO
  ENDIF

END SUBROUTINE init_fft_storage_scheme
