SUBROUTINE peak_div_u(u,peak)
  USE defprecision_module
  USE state_module, ONLY:velocity,ltime0
  USE diff_op_module, ONLY:div
  USE mpi_transf_module, ONLY:mysx_spec,myex_spec,mysy_spec,myey_spec
  USE parameter_module, ONLY: Nmax  
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  TYPE(velocity) :: u
  COMPLEX(kind=kr),POINTER :: divu(:,:,:)
  REAL(kind=kr) :: peak,lpeak
  INTEGER :: ierr

  ALLOCATE(divu(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
  divu = div(u%spec(:,:,:,:,ltime0))
!  divu = div(u%rhs (:,:,:,:))
  lpeak = MAXVAL(ABS(divu))
!  print*,maxloc(abs(divu))
  DEALLOCATE(divu)
  CALL MPI_REDUCE(lpeak,peak,1,PM_MPI_FLOAT_TYPE,MPI_MAX, &
       &          0,MPI_COMM_WORLD,ierr)

END SUBROUTINE peak_div_u
