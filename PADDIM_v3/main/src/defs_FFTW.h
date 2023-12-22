!/************************************************************/
!/**  THIS INCLUDE FILE DEFINES WHICH FFTW SUBROUTINES ARE  **/
!/**      CALLED FOR SINGLE AND DOUBLE PRECISION            **/
!/************************************************************/

#if defined(SINGLE_PRECISION)
#define PM_FFTW_EXECUTE_DFT sfftw_execute_dft
#define PM_FFTW_EXECUTE_DFT_C2R sfftw_execute_dft_c2r
#define PM_FFTW_EXECUTE_DFT_R2C sfftw_execute_dft_r2c
#define PM_FFTW_PLAN_MANY_DFT_R2C sfftw_plan_many_dft_r2c
#define PM_FFTW_PLAN_MANY_DFT_C2R sfftw_plan_many_dft_c2r
#define PM_FFTW_PLAN_MANY_DFT sfftw_plan_many_dft
#define PM_FFTW_CLEANUP sfftw_cleanup
#endif
#if defined(DOUBLE_PRECISION) 
#define PM_FFTW_EXECUTE_DFT dfftw_execute_dft
#define PM_FFTW_EXECUTE_DFT_C2R dfftw_execute_dft_c2r
#define PM_FFTW_EXECUTE_DFT_R2C dfftw_execute_dft_r2c
#define PM_FFTW_PLAN_MANY_DFT_R2C dfftw_plan_many_dft_r2c
#define PM_FFTW_PLAN_MANY_DFT_C2R dfftw_plan_many_dft_c2r
#define PM_FFTW_PLAN_MANY_DFT dfftw_plan_many_dft
#define PM_FFTW_CLEANUP dfftw_cleanup
#endif
