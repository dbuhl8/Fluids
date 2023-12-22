! This subroutine computes the dissipation rate of buoyancy fields.
! For a scalar field x, given in spectral space, the routine computes 
! < ( \nabla x) ** 2 > , where <...> denotes a spatial average.  
! Not implemented efficiently, needs further optimization if called frequently.
FUNCTION dissipation_buo(x)
  USE defprecision_module
  USE mpi_transf_module, ONLY:mysx_spec,myex_spec,mysy_spec,myey_spec
  USE parameter_module, ONLY: Nmax
  USE diff_op_module, ONLY: d_by_dx,d_by_dy,d_by_dz
  IMPLICIT NONE
  REAL(kind=kr)    :: dissipation_buo
  COMPLEX(kind=kr) :: x(:,:,:)
  COMPLEX(kind=kr) :: work(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec)

  dissipation_buo = 0._kr
  
  work = d_by_dx(x)
  dissipation_buo = (rms(work))**2

#ifndef TWO_DIMENSIONAL
  work = d_by_dy(x)
  dissipation_buo = dissipation_buo + (rms(work))**2
#endif

  work = d_by_dz(x)
  dissipation_buo = dissipation_buo + (rms(work))**2
  
END FUNCTION dissipation_buo
