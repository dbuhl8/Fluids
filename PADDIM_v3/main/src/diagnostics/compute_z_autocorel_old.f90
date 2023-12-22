subroutine compute_z_autocorel(vel,autocorel)
use defprecision_module
use defs_2D_3D_module
  USE message_passing_module, ONLY : myid
  use parameter_module, ONLY : Nz,Gammax,Gammay,Gammaz,Nmax,Lmax,kz
  USE mpi_transf_module, ONLY:  mysy_spec,myey_spec,mysx_spec,myex_spec
  USE MPI
  implicit none

  real(kind=kr), allocatable :: lauto_spec(:),auto_spec(:,:)
  complex (kind=kr) :: vel(0:,mysx_spec:,mysy_spec:,vec_x:)
  real(kind=kr) :: autocorel(0:Nz-1,vec_x:vec_z)
  real(kind=kr) :: dz,zc,fac


  integer :: i,j,k,n,ierr

! Returns autocorrelation function in the root proc 


  allocate(lauto_spec(0:2*Nmax-1))
  allocate(auto_spec(0:2*Nmax-1,vec_x:vec_z))
   
  auto_spec = 0.

! Compute local power summed over kx,ky
  do n = vec_x,vec_z
     do k=0,2*Nmax-1   ! This contains everything from -N_z/2+1 to Nz/2
        lauto_spec(k)= 0.
        do j=mysy_spec,myey_spec
           do i=mysx_spec,myex_spec
             fac = REAL((2 - DIM(1,i)),kr) 
!           if(i.eq.0.or.i.eq.Lmax) then
!                   fac = 1.
!           else
!                   fac = 2.
!           endif
             lauto_spec(k)=lauto_spec(k)+vel(k,i,j,n)*CONJG(vel(k,i,j,n))*fac 
         enddo
       enddo
     enddo

! Add them all into the auto_spec of root comm.
     call MPI_reduce(lauto_spec,auto_spec(0,n),2*Nmax,PM_MPI_FLOAT_TYPE,MPI_SUM,0,  &
          &   MPI_COMM_WORLD,ierr)
  enddo

  autocorel = 0.

  if(myid.eq.0) then

     dz = Gammaz/Nz

     do n=vec_x,vec_z
        do i = 0,Nz/2
          autocorel(i,n) = auto_spec(0,n)
          zc = i*dz
          do k=1,2*Nmax-1
              autocorel(i,n) = autocorel(i,n)+cos(kz(k)*zc)*auto_spec(k,n)
          enddo
        enddo
     enddo
  endif

  deallocate(auto_spec,lauto_spec)

end subroutine compute_z_autocorel

