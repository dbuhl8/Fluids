! This is the function used for error handling.
SUBROUTINE pn_check(status)
  USE defprecision_module
  IMPLICIT NONE
#ifdef BLUE_GENE
  INCLUDE 'pnetcdf.inc'
#else
#include "pnetcdf.inc"
#endif
  INTEGER, INTENT ( in) :: status
    
  IF(status /= nf_noerr) THEN
     PRINT *,status
     STOP "Stopped"
  END IF
END SUBROUTINE pn_check
