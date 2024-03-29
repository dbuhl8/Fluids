subroutine compute_power_input(power,u)
  use defprecision_module
  use parameter_module, ONLY : Nx,Ny,Nz,Gammax,Gammay,Gammaz
  use state_module, ONLY: velocity
  USE mpi_transf_module, ONLY:  mysy_phys,myey_phys,mysz_phys,myez_phys
  use message_passing_module, ONlY: myid
  use forcing_module, ONLY: str, force_real
  use MPI
  implicit none
  real(kind=kr)  :: power
  type(velocity) :: u
  real (kind=kr) :: lpower, dx,dy,dz,dv,volume
  integer(kind=ki) :: i,j,k,l,ierr
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

!  print *, "cpu "//trim(str(myid))//": starting my power loop"
  lpower = 0._kr
  do l = vec_x, vec_z
     do k=mysz_phys,myez_phys
        do j=mysy_phys,myey_phys
            do i=0,Nx-1 ! pencils are complete along x dir in phys space, and complete along z in spec space
              lpower = lpower + dv*u%phys(i,j,k,l)*force_real(i,j,k,l)
           enddo
        enddo
     enddo
  enddo
  !print *, "cpu "//trim(str(myid))//": got past lpower loop"

  lpower = lpower / volume 

  call MPI_REDUCE(lpower,power,1,PM_MPI_FLOAT_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  !print *, "cpu "//trim(str(myid))//": got past mpi reduce"
  
end subroutine compute_power_input
