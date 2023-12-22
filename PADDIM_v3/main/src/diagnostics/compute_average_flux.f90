subroutine compute_average_flux(flux,buo,u)
  use defprecision_module
  use parameter_module, ONLY : Nx,Ny,Nz,Gammax,Gammay,Gammaz
  use state_module, ONLY: buoyancy,velocity
  USE mpi_transf_module, ONLY:  mysy_phys,myey_phys,mysz_phys,myez_phys
  use MPI
  implicit none
  real(kind=kr)  :: flux
  type(buoyancy) :: buo
  type(velocity) :: u
  real (kind=kr) :: lflux,dx,dy,dz,dv,volume
  integer(kind=ki) :: i,j,k,ierr
  dx = Gammax / Nx
  dz = Gammaz / Nz
#ifdef TWO_DIMENSIONAL
  dv = dx*dz
  volume = Gammax * Gammaz
#else
  dy = Gammay / Ny
  dv = dx*dy*dz
  volume = Gammax * Gammay * Gammaz
#endif

  lflux = 0._kr
  do k=mysz_phys,myez_phys
     do j=mysy_phys,myey_phys
        do i=0,Nx-1
           lflux = lflux + dv*u%phys(i,j,k,vec_z)*buo%phys(i,j,k)
        enddo
     enddo
  enddo

  lflux = -lflux / volume 

  call MPI_REDUCE(lflux,flux,1,PM_MPI_FLOAT_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  
end subroutine compute_average_flux
