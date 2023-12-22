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

FUNCTION rms_spec_scalar(in)
  USE defprecision_module
  USE parameter_module, ONLY  :  Nmax
  USE mpi_transf_module,  ONLY:  mysx_spec,myex_spec,mysy_spec,myey_spec
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  REAL (kind=kr)              :: rms_spec_scalar
  COMPLEX (kind=kr)           :: in(0:,mysx_spec:,mysy_spec:)
  INTEGER (kind=ki)           :: k,i,j,ierr
  REAL(kind=kr)               :: lenergy,genergy,fac
!
  lenergy=0._kr
  DO j=mysy_spec,myey_spec
     DO i=mysx_spec,myex_spec
        fac = REAL((2 - DIM(1,i)),kr) 
        DO k=0,2*Nmax-1
           lenergy = lenergy + in(k,i,j) * CONJG(in(k,i,j)) * fac 
        ENDDO
     ENDDO
  ENDDO
  genergy = 0._kr
  CALL MPI_REDUCE(lenergy,genergy,1,PM_MPI_FLOAT_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  rms_spec_scalar = SQRT(genergy)
!
END FUNCTION rms_spec_scalar
