SUBROUTINE init_u_phys(u)
  USE defprecision_module
  USE parameter_module, ONLY: &
     Nx,Ny,Nz,Gammax,Gammay,Gammaz,pi
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys
  IMPLICIT NONE
  TYPE(velocity) :: u
  INTEGER(kind=ki) :: i,j,k,idum
  REAL(kind=kr) :: xc,yc,zc,dx,dy,dz,rn

  dx = Gammax / Nx
  dy = Gammay / Ny
  dz = Gammaz / Nz 
  idum = -7

  DO k=mysz_phys,myez_phys
     zc = k*dz
     DO j=mysy_phys,myey_phys
        yc = j*dy
        DO i=0,Nx-1
           xc = i*dx
#ifdef TWO_DIMENSIONAL
           ! Initial velocity for 2D case 
           u%phys(i,j,k,vec_x) = 0       !sin(yc)    !1E-15*rn
           rn = decomp_independent_random(i,j,k,idum) - 0.5_kr
           u%phys(i,j,k,vec_z) = 0.    !1E-15*rn
#else
           ! Initial velocity for 3D case 
           !rn = decomp_independent_random(i,j,k,idum) - 0.5_kr
           u%phys(i,j,k,vec_x) = 0.       !sin(yc)
           rn = decomp_independent_random(i,j,k,idum) - 0.5_kr
           u%phys(i,j,k,vec_y) =  0.
           u%phys(i,j,k,vec_z) = 1.0d-4*rn
#endif
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE init_u_phys
