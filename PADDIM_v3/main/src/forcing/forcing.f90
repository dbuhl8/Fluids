!File: forcing.f90
!Author: Dante Buhl, Pascale Garaud
!Purpose: return a spectral and physical forcing array to crhs_velocity

SUBROUTINE forcing(t)
  ! This routine takes as an input, time, and returns forcing in spectral and physical forms
  USE defprecision_module
  USE parameter_module
  USE message_passing_module, ONLY: myid
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys,myex_phys, mysx_phys, &
      &                        mysx_spec,myex_spec,mysy_spec,myey_spec, mysz_spec, myez_spec, &
      &                        FFT_r2c, FFT_c2r
  IMPLICIT NONE
  REAL (kind=kr) :: hkx,hky,hkz
  REAL (kind=kr) :: ksquare
  REAL (kind=kr), dimension(2) :: gp
  REAL (kind=kr) :: xfactor, yfactor, scaling
  INTEGER (kind=ki) :: i,j,k,l
  REAL(kind=kr) :: yc,dy
  REAL(kind=kr) :: t

! This routine checks whether we want to define the forcing in real or complex space and then cast it into the other
  ! this can be done in a different file 

  if (c2r) then
    force_spec = (0._kr, 0._kr)
    force_real = 0._kr
    if (forcedCPU) then
      ! This implements a STOCHASTIC_FORCING method 
      scaling = 7./45.
      do j=mysy_spec,myey_spec
        hky = ky(j)
        do i=mysx_spec,myex_spec
          hkx = kx(i)
          ksquare=sqrt(hkx**2 + hky**2)
          if((ksquare .le. kmax_forcing) .and. (hkx .ne. 0 .and. hky .ne. 0)) then
            ksquare = max(ksquare,epsilon(1._kr)) ! avoid floating exception in 1/ksquare later...
            gp = gpinterp(t, i, j)
            xfactor = hky*scaling/ksquare
            yfactor = -hkx*scaling/ksquare
            force_spec(0, i, j, vec_x) = xfactor*gp(1)*(1._kr, 0._kr) + xfactor*gp(2)*(0._kr, 1._kr)
            force_spec(0, i, j, vec_y) = yfactor*gp(1)*(1._kr, 0._kr) + yfactor*gp(2)*(0._kr, 1._kr)
          end if
        enddo
      enddo
    end if
   
    ! Spectral forcing is cast into real forcing
    do i = vec_x, vec_z
      CALL FFT_c2r(force_spec(:,:,:,i), force_real(:,:,:,i))
    enddo
  else
    ! real space forcing
    dy = gammay / ny
    do k=mysz_phys, myez_phys            
      do j=mysy_phys,myey_phys
        do i=mysx_phys, myex_phys
          yc = j*dy
          !place the forcing here
          force_real(i, j, k, vec_x) = sin(yc)  !a sinusoidal shear flow
          force_real(i, j, k, vec_y) = 0._kr   
          force_real(i, j, k, vec_z) = 0._kr
        enddo
      enddo
    enddo
                      
    ! real forcing is cast into spectral forcing                                                           
    do i=vec_x, vec_z
        call fft_r2c(force_real(:,:,:,i), force_spec(:,:,:,i))
    enddo
  end if

  
END SUBROUTINE forcing

