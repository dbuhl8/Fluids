!+--------------------------------------------------------------------+
!| The following function computes the volume integral over a 3d      |
!| scalar field:                                                      |
!|                                                                    |
!| volume_integral(\phi) = \int\int\int phi(x,y,z) dx dy dz           |
!|                                                                    |
!| where the integral is taken over the domain                        |
!| [0,Gammax] \times [0,Gammay] \times [0,Gammaz]                     |
!| and \phi denotes the scalar field.                                 |
!| The result is returned on process 0 in MPI_COMM_WORLD only.        |
!+--------------------------------------------------------------------+
FUNCTION volume_integral(in)
  USE defprecision_module
  USE parameter_module, ONLY  : Nx,Ny,Nz,Gammax,Gammay,Gammaz
  USE mpi_transf_module, ONLY:  mysy_phys,mysz_phys,mysz_phys,myez_phys
  USE message_passing_module, ONLY : myid
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  REAL (kind=kr)     :: volume_integral
  REAL (kind=kr)     :: in(0:,mysy_phys:,mysz_phys:)
  REAL (kind=kr)     :: local_int,dx,dy,dz,dv
  INTEGER(kind=ki)   :: ierr

  ! Use trapezoidal rule to evaluate integrals.
  ! For periodic functions, this has spectral accuracy (see Boyd, 2000, §19.8.4).

  
  dx = Gammax / Nx
  dz = Gammaz / Nz
#ifdef TWO_DIMENSIONAL
  dv = dx*dz
#else
  dy = Gammay / Ny
  dv = dx*dy*dz
#endif

  local_int = dv*SUM(in)
  
  CALL MPI_REDUCE(local_int,volume_integral,1,PM_MPI_FLOAT_TYPE,MPI_SUM, &
       &          0,MPI_COMM_WORLD,ierr)

END FUNCTION volume_integral
