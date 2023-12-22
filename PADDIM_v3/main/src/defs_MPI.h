!/************************************************************/
!/**  THIS INCLUDE FILE DEFINES WHICH MPI DATA TYPES ARE    **/
!/**      USED FOR FLOATING POINT COMPUTATIONS              **/
!/************************************************************/

#if defined(SINGLE_PRECISION)
#define PM_MPI_FLOAT_TYPE MPI_REAL    
#define PM_MPI_COMPLEX_TYPE MPI_COMPLEX    
#endif
#if defined(DOUBLE_PRECISION) 
#define PM_MPI_FLOAT_TYPE MPI_DOUBLE_PRECISION    
#define PM_MPI_COMPLEX_TYPE MPI_DOUBLE_COMPLEX    
#endif
#define PM_MPI_SINGLE_FLOAT_TYPE MPI_REAL 
#define PM_MPI_DOUBLE_FLOAT_TYPE MPI_DOUBLE_PRECISION 
