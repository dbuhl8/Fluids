subroutine error_Temp_adv_diff(err,Temp,t)
  use defprecision_module
  use state_module, ONLY : buoyancy
  use mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys
  use parameter_module, ONLY: Nx,Ny,Nz,Gammax,Gammay,Gammaz,alpha,beta,gamma
  use diagnostics_module, ONLY: rms
  implicit none
  type(buoyancy) :: Temp
  real(kind=kr) :: err,t
  Real(kind=kr) :: Tempa(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys)
  integer(kind=ki) :: i,j,k
  real(kind=kr) :: xc,yc,zc,dx,dy,dz,decay_fac

  dx = Gammax / Nx
  dy = Gammay / Ny
  dz = Gammaz / Nz 

  decay_fac= exp( -(alpha**2 + (3*beta)**2 + (2*gamma)**2)*t )
  do k=mysz_phys,myez_phys
     zc = k*dz
     do j=mysy_phys,myey_phys
        yc = j*dy
        do i=0,Nx-1
           xc = i*dx
           Tempa(i,j,k) = decay_fac*cos(alpha*(xc-100._kr*t)+3*beta*yc+2*gamma*zc)
        enddo
     enddo
  enddo

  err = rms(Temp%phys - Tempa)
  
end subroutine error_Temp_adv_diff
