! This subroutine computes the dissipation rate of buoyancy fields.
! For a scalar field x, given in spectral space, the routine computes 
! < ( \nabla x) ** 2 > , where <...> denotes a spatial average.  
! No implemented efficiently, needs further optimization if called frequently.
function dissipation_buo(x)
  use defprecision_module
  use mpi_transf_module, ONLY:mysx_spec,myex_spec,mysy_spec,myey_spec
  use diff_op_module, ONLY: d_by_dx,d_by_dy,d_by_dz
  implicit none
  real(kind=kr)    :: dissipation_buo
  complex(kind=kr) :: x(:,:,:)
  complex(kind=kr) :: work(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec)

  dissipation_buo = 0._kr
  
  work = d_by_dx(x)
  dissipation_buo = (rms(work))**2

  work = d_by_dy(x)
  dissipation_buo = dissipation_buo + (rms(work))**2

  work = d_by_dz(x)
  dissipation_buo = dissipation_buo + (rms(work))**2
  
end subroutine compute_cuo_dissipation
