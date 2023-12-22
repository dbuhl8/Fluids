! This subroutine closes the netCDF file used to store simulation data
SUBROUTINE pn_close_simdat_file
  USE defprecision_module
  IMPLICIT NONE
#include "pnetcdf.inc"
  CALL pn_check( nfmpi_close(ncid_simdat) )
END SUBROUTINE pn_close_simdat_file
