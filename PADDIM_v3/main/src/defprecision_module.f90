MODULE defprecision_module
! Definition of the precision of the data types used 
! This is compiler specific !!!
! The precision of the datatypes has to correspont to 
! the precision expected by the FFT and MPI libraries used.
! The code can be compiled in single or double precision accuracy - see cpp directives 
! below. The include files defs_FFTW.h and defs_MPI.h ensure that the right routines 
! and MPI datatypes are used. The definitions given in these files and in this 
! module must therefore be consistent!!!
  IMPLICIT NONE
  SAVE
#ifdef SINGLE_PRECISION
  INTEGER, PARAMETER :: kr      = kind(0.)   !single precision real kind 
#endif
#ifdef DOUBLE_PRECISION
  INTEGER, PARAMETER :: kr      = kind(dble(0.))   !double precision real kind 
#endif
  INTEGER, PARAMETER :: krs     = kind(0.)   !single precision real kind, used for example
                                             !for writing some output files
  INTEGER, PARAMETER :: krd     = kind(dble(0.)) !double precision real kind, may be used for 
                                                 !code regions prone to round-off errors
  INTEGER, PARAMETER :: ki      = kind(1)    !integer kind used
  INTEGER, PARAMETER :: kli     = selected_int_kind(18) !long integer kind used
                                  ! Should be at least Nx*Ny*Nz, so choose an 
                                  ! 8 byte integer  
  INTEGER, PARAMETER :: kiMPI   = ki    !integer kind used for MPI constants
  INTEGER, PARAMETER :: kr_jc   = krs   !real kind used for JC-Files 
 
END MODULE defprecision_module
