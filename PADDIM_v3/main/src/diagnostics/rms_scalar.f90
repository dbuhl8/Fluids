!+--------------------------------------------------------------------+
!| The following function computes the L2-norm of a 3d scalar field.  |
!| If the field is denotes by \phi, the L2-norm is defined as         |
!|                                                                    |
!| rms_scalar(\phi) =                                                 |
!|            sqrt( 1/V * \int\int\int \phi(x,y,z)**2 dx dy dz )      |
!|                                                                    |
!| where the integral is taken over the domain                        |
!| [0,Gammax] \times [0,Gammay] \times [0,Gammaz]                     |
!| and V = Gammax * Gammay * Gammaz.                                  |
!| IMPORTANT NOTE: The implementation is readable, but probably quite |
!| inefficient. Will need further optimization later...               |
!+--------------------------------------------------------------------+
FUNCTION rms_scalar(in)
  USE defprecision_module
  USE parameter_module, ONLY  :  Nx,Ny,Nz,Gammax,Gammay,Gammaz
  USE integral_module, ONLY : volume_integral
  USE mpi_transf_module,  ONLY:  mysy_phys,myey_phys,mysz_phys,myez_phys
  IMPLICIT NONE
  REAL (kind=kr)              :: rms_scalar
  REAL (kind=kr)              :: in(:,:,:)
  REAL (kind=kr), POINTER     :: work(:,:,:)
  REAL (kind=kr)              :: volume
!
#ifdef TWO_DIMENSIONAL
  volume=Gammaz*Gammax
#else 
  volume=Gammaz*Gammax*Gammay
#endif
  ALLOCATE(work(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys))
  work = in**2
  rms_scalar=SQRT( volume_integral(work)/volume )
  DEALLOCATE(work)
!
END FUNCTION rms_scalar
