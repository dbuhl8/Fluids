!+--------------------------------------------------------------------+
!| This subroutine performs initializations needed to compute         |
!| transforms from physical space to fourier space and back.          |
!|                                                                    |
!| INPUT:                                                             |
!|                                                                    |
!|   maxL                 = highest frequency of x-modes              |
!|   maxM                 = highest frequency of y-modes              |
!|   maxN                 = highest frequency of z-modes              |
!|   numx                 = number of grid points in x-direction      |
!|                          ( numx >= 2*maxL, must be even )          |
!|   numy                 = number of grid points in y-direction      |
!|                          ( numx >= 2*maxM, must be even )          |
!|   numz                 = number of grid points in z-direction      |
!|                          ( numx >= maxN, must be even )            |
!|   nprocs1              = number of processes for y-decomposition   |
!|                          in physical space (pencil decomposiotion) |
!|   nprocs2              = number of processes for z-decomposition   |
!|                          in physical space (pencil decomposiotion) |
!|   root_communicator    = MPI communicator to be used               |
!|   fftw_wisdom_filename = filename of the file in which FFTW should |
!|                          store wisdom information. If omitted,     |
!|                          no FFTW wisdom is stored or read from     |
!|                          a disk file.                              |
!|   fftw_wisdom_fileunit = unit to which the FFTW wisdom file should |
!|                          be connected                              |
!|                                                                    |
!| If numx > 2*maxL or numy > 2*maxM or numz > 2*maxN, the spectra    |
!| are  padded with zeros before a transform from spectral to         |
!| physical space is performed.                                       |
!| For a real to complex transfrom, only coefficients for             |
!| (l,m,n) with 0 <= l <=Lmax, -Mmax <= m <= Mmax, -Nmax <= n <=Nmax  |
!| are returned. For l=Lmax, |m|=Mmax and |n|=Nmax, the output is     |
!| undefined. These modes are typically not used in pseudospectral    |
!| codes.                                                             |
!|                                                                    |
!+--------------------------------------------------------------------+
!| Author: Stephan Stellmach                                          |
!+--------------------------------------------------------------------+
SUBROUTINE init_transforms(maxL,maxM,maxN,                          &
&                          numx,numy,numz,                          &
&                          nprocs1,nprocs2,root_communicator,       &
&                          fftw_wisdom_filename,fftw_wisdom_fileunit)
  USE defprecision_module
  IMPLICIT NONE
  INTEGER(kind=ki)            :: maxL,maxM,maxN
  INTEGER(kind=ki)            :: numx,numy,numz
  INTEGER(kind=ki)            :: nprocs1,nprocs2
  INTEGER(kind=kiMPI)         :: root_communicator
  CHARACTER(LEN=100),OPTIONAL :: fftw_wisdom_filename 
  INTEGER(kind=ki),OPTIONAL   :: fftw_wisdom_fileunit 
  INTEGER(kind=kiMPI)         :: myid_root_communicator,ierr

  ! store size of global transform
  Lmax = maxL
  Mmax = maxM
  Nmax = maxN 

  Nx=numx
  Ny=numy
  Nz=numz

#ifdef TWO_DIMENSIONAL
  ! For 2d case, check if numy=1 and maxM=0
  IF (maxM.NE.0 .OR. Numy.NE.1) THEN 
     PRINT *,"Error: Ny or Mmax not compatible with 2d version"
     STOP
  ENDIF
#endif

  ! check if numx,numy,numz are all even - except in the 2d case, where numy = 1 is allowed
  IF (MOD(numx,2).NE.0 .OR. (MOD(numy,2).NE.0 .AND. numy.NE.1) .OR. MOD(numz,2).NE.0) THEN 
     PRINT*,"Error: Number of grid points in at least one direction not even..."
     CALL MPI_FINALIZE(ierr)
     STOP
  ENDIF

  ! check if numx >= 2*Lmax, numy >= 2*Mmax, numz >= 2*Nmax
  IF (numx.LT.2*Lmax .OR. numy.LT.2*Mmax .OR. numz.LT.2*Nmax) THEN 
     PRINT*,"Error: Number of grid points to low for given spectral resolution..."
     CALL MPI_FINALIZE(ierr)
     STOP
  ENDIF

  ! try to import information stored after previsous initializations of FFTW
  IF (PRESENT(fftw_wisdom_filename)) THEN 
     IF (.NOT.PRESENT(fftw_wisdom_fileunit)) &
     &  PRINT*,"Error: an fftw wisdom file unit must be specified!"
     use_FFTW_wisdom_file = .TRUE.
     FFTW_wisdom_file     = fftw_wisdom_filename
     uwisdom              = fftw_wisdom_fileunit
     CALL MPI_COMM_RANK(root_communicator,myid_root_communicator,ierr)
     WRITE(FFTW_wisdom_file(LEN_TRIM(FFTW_wisdom_file)+1:),'(a,I5.5)') &
          &   "_",myid_root_communicator
  ELSE 
     use_FFTW_wisdom_file = .FALSE.
  ENDIF
  CALL get_FFTW_wisdom

#ifdef TWO_DIMENSIONAL
  CALL init_transf_2d(root_communicator)
#else
  CALL init_transf_pencil(root_communicator,nprocs1,nprocs2)
#endif

! save FFTW wisdom to file
  CALL save_FFTW_wisdom

END SUBROUTINE init_transforms

