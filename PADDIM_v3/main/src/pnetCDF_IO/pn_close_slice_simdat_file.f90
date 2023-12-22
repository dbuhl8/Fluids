! This subroutine closes the 3 netCDF file used to store slice data
SUBROUTINE pn_close_slice_simdat_file
  USE defprecision_module
  IMPLICIT NONE
#include "pnetcdf.inc"
  CALL pn_check( nfmpi_close(ncid_xyslice_simdat) )
  CALL pn_check( nfmpi_close(ncid_xzslice_simdat) )
  CALL pn_check( nfmpi_close(ncid_yzslice_simdat) )

END SUBROUTINE pn_close_slice_simdat_file
