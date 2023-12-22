!+--------------------------------------------------------------------+
!| The following function computes the L2-norm of a 3d scalar field.  |
!| The field is assumed to be given in spectral space.                |
!| If the field is denotes by \phi, the L2-norm is defined as         |
!|                                                                    |
!| rms_scalar(\phi) =                                                 |
!|            sqrt( 1/V * \int\int\int \phi(x,y,z)**2 dx dy dz )      |
!|                                                                    |
!| where the integral is taken over the domain                        |
!| [0,Gammax] \times [0,Gammay] \times [0,Gammaz]                     |
!| and V = Gammax * Gammay * Gammaz.                                  |
!+--------------------------------------------------------------------+

FUNCTION rms_spec_vector(in)
  USE defprecision_module
  USE parameter_module, ONLY  :  Nmax
  USE integral_module, ONLY : volume_integral
  USE mpi_transf_module,  ONLY:  mysx_spec,myex_spec,mysy_spec,myey_spec
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  REAL (kind=kr)              :: rms_spec_vector
  complex (kind=kr)           :: in(0:,mysx_spec:,mysy_spec:,vec_x:)
  integer (kind=ki)           :: k,i,j,n,ierr
  real(kind=kr)               :: lenergy,genergy,fac
!
  lenergy=0._kr

  do n=vec_x,vec_z
     do j=mysy_spec,myey_spec
        do i=mysx_spec,myex_spec
           fac = REAL((2 - DIM(1,i)),kr) 
           do k=0,2*Nmax-1
              lenergy = lenergy + in(k,i,j,n) * CONJG(in(k,i,j,n)) * fac 
           enddo
        enddo
     enddo
  enddo
  genergy = 0._kr
  call MPI_REDUCE(lenergy,genergy,1,PM_MPI_FLOAT_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  rms_spec_vector = sqrt(genergy)
!
END FUNCTION rms_spec_vector
