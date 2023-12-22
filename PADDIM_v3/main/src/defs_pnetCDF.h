!/****************************************************************/
!/**  THIS INCLUDE FILE DEFINES WHICH PNETCDF SUBROUTINES ARE   **/
!/**          CALLED FOR SINGLE AND DOUBLE PRECISION            **/
!/****************************************************************/

#if defined(SINGLE_PRECISION)
#define PM_NF_FLOAT NF_real
#define PM_NFMPI_PUT_VAR_FLOAT nfmpi_put_var_real
#define PM_NFMPI_GET_VAR_FLOAT nfmpi_get_var_real
#define PM_NFMPI_GET_VAR_FLOAT_ALL nfmpi_get_var_real_all
#define PM_NFMPI_PUT_VARA_FLOAT nfmpi_put_vara_real
#define PM_NFMPI_PUT_VARA_FLOAT_ALL nfmpi_put_vara_real_all
#define PM_NFMPI_GET_VARA_FLOAT_ALL nfmpi_get_vara_real_all

#endif
#if defined(DOUBLE_PRECISION) 
#define PM_NF_FLOAT NF_double
#define PM_NFMPI_PUT_VAR_FLOAT nfmpi_put_var_double
#define PM_NFMPI_GET_VAR_FLOAT nfmpi_get_var_double
#define PM_NFMPI_GET_VAR_FLOAT_ALL nfmpi_get_var_double_all
#define PM_NFMPI_PUT_VARA_FLOAT nfmpi_put_vara_double
#define PM_NFMPI_PUT_VARA_FLOAT_ALL nfmpi_put_vara_double_all
#define PM_NFMPI_GET_VARA_FLOAT_ALL nfmpi_get_vara_double_all
#endif
