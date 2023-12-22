!File: forcing.f90
!Author: Dante Buhl, Pascale Garaud
!Purpose: return a spectral forcing array for crhs_velocity



SUBROUTINE forcing(t, force_spec)
! This routine takes as an input, time, and returns forcing in spectral space
  USE defprecision_module
  USE parameter_module
  USE message_passing_module, ONLY: myid
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys,myex_phys, mysx_phys, &
      &                        mysx_spec,myex_spec,mysy_spec,myey_spec, mysz_spec, myez_spec, &
      &                        FFT_r2c
  IMPLICIT NONE
  COMPLEX(kind=kr) :: force_spec(0:2*Nmax-1, mysx_spec:myex_spec, mysy_spec:myey_spec, 1:3) ! RHS
  REAL (kind=kr),POINTER     :: force_real(:,:,:,:)
  REAL (kind=kr) :: hkx,hky,hkz
  REAL (kind=kr) :: ksquare, gpholder
  REAL (kind=kr) :: xfactor, yfactor
  INTEGER (kind=ki) :: i,j,k,l
  REAL(kind=kr) :: yc,dy
  REAL(kind=kr) :: t

 
#ifdef STOCHASTIC_FORCING

    print *, "cpu "//trim(str(myid))//": code is about to get through forcing loop"
    force_spec = (0._kr, 0._kr)
    
    if(myid .eq. 112) then
        print *, "cpu "//trim(str(myid))//":  mysx:", mysx_spec, ", myex:", myex_spec, ", mysy:", mysy_spec, &
        & " myey_spec:", myey_spec, ", mysz:", mysz_spec, " myez_spec:", myez_spec, ", 2*nmax-1:, ", 2*Nmax-1
    end if
    DO j=mysy_spec,myey_spec
        hky = ky(j)
        DO i=mysx_spec,myex_spec
            hkx = kx(i)
            ksquare=sqrt(hkx**2 + hky**2)
            if((ksquare .le. KMAX_forcing) .and. (ksquare .ne. 0)) then
                ksquare = MAX(ksquare,EPSILON(1._kr)) ! Avoid floating exception in 1/ksquare later...
                if(myid .eq. 112) then
                    print *, "cpu "//trim(str(myid))//":  i", i, ", j", j
                end if
                gpholder = gpinterp(t, i, j)
                xfactor = hky*gpholder/ksquare
                yfactor = -hkx*gpholder/ksquare
                force_spec(0, i, j, vec_x) = xfactor*(1._kr, 0._kr) 
                force_spec(0, i, j, vec_y) = yfactor*(1._kr, 0._kr) 
            end if
        ENDDO
    ENDDO
    print *, "cpu "//trim(str(myid))//": code went through forcing loop"

#else

    ALLOCATE(force_real(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys,vec_x:vec_z))
    dy = Gammay / Ny
    DO k=mysz_phys, myez_phys            
        DO j=mysy_phys,myey_phys
            DO i=mysx_phys, myex_phys
                yc = j*dy
                force_real(i, j, k, vec_x) = sin(yc) 
                force_real(i, j, k, vec_y) = 0._kr   
                force_real(i, j, k, vec_z) = 0._kr
            ENDDO
        ENDDO
    ENDDO
                                                                                 
    DO i=vec_x, vec_z
        CALL FFT_r2c(force_real(:,:,:,i), force_spec(:,:,:,i))
    ENDDO  
    DEALLOCATE(force_real)

#endif
    
    contains

    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

END SUBROUTINE forcing

