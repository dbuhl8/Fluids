!+--------------------------------------------------------------------+
!| The following function computes the L2-norm of a 3d vector field.  |
!| If the field is denoted by \phi, the L2-norm is defined as         |
!|                                                                    |
!| rms_scalar(\phi) =                                                 |
!| sqrt( 1/V * \int\int\int \phi_x^2 + \phi_y^2 + \phi_z^2 dx dy dz ) |
!|                                                                    |
!| where the integral is taken over the domain                        |
!| [0,Gammax] \times [0,Gammay] \times [0,Gammaz]                     |
!| and V = Gammax * Gammay * Gammaz.                                  |
!| NOTE: The implementation here is very readable, but likely to be   |
!| inefficient. Will need further optimization....                    |  
!+--------------------------------------------------------------------+
FUNCTION rms_vector(in)
  USE defprecision_module
  USE parameter_module, ONLY  :  Nx,Ny,Nz,Gammax,Gammay,Gammaz
  USE integral_module, ONLY : volume_integral
  USE mpi_transf_module,  ONLY:  mysy_phys,myey_phys,mysz_phys,myez_phys
  IMPLICIT NONE
  REAL (kind=kr)              :: rms_vector
  REAL (kind=kr)              :: in(:,:,:,:)
  REAL (kind=kr),ALLOCATABLE  :: work(:,:,:)
  REAL (kind=kr)              :: volume
!
  ALLOCATE(work(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys))
#ifdef TWO_DIMENSIONAL
  volume=Gammaz*Gammax
  work = in(:,:,:,vec_x)**2 + in(:,:,:,vec_z)**2
#else
  volume=Gammaz*Gammax*Gammay
  work = in(:,:,:,vec_x)**2 + in(:,:,:,vec_y)**2 + in(:,:,:,vec_z)**2
#endif
  rms_vector=SQRT( volume_integral(work)/volume )
  DEALLOCATE(work)
!
END FUNCTION rms_vector
