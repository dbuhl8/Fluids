SUBROUTINE pn_close_dump
  USE defprecision_module
  IMPLICIT NONE
#include "pnetcdf.inc" 

  CALL pn_check( nfmpi_close(ncid_dump) )

END SUBROUTINE pn_close_dump
