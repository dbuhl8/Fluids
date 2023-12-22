! The following function assumes that a global field is initialized with 
! random numbers in the following way:
!
! idummy = idum
! do k=0,Nz
!   do j=0,Ny-1
!     do i=0,Nx-1
!       x(i,j,k) = ran2(idummy)
!     enddo
!   enddo
! enddo
!
! where idum is the initialization integer of the routine ran2 from the
! Numerical Recipes in F77 book. 
! The following function can be used to compute x(i,j,k).
! This allows to initialize a field with random noise in a portable, reproducible
! and decomposition independent way. idum has to be smaller than zero.
FUNCTION decomp_independent_random(i,j,k,idum)
  USE defprecision_module
  USE parameter_module, ONLY : Nx,Ny
  IMPLICIT NONE
  REAL(kind=kr) :: decomp_independent_random
  INTEGER(kind=ki) :: i,j,k
  INTEGER(kind=ki) :: num_of_evals_total = 0
  INTEGER(kind=ki) :: idum 
  INTEGER(kind=ki) :: idum_init = 0
  INTEGER(kind=ki) :: jdum = 0
  INTEGER(kind=ki) :: num_of_evals,count
  REAL(kind=kr) :: rbla !,ran2

  ! check if things should be initialized again 
  ! (i.e. idum changed since last call)
  IF (idum .NE. idum_init) THEN 
     jdum = idum 
     idum_init = idum
     num_of_evals_total = 0 
  ENDIF

  num_of_evals = (k*(Nx*Ny) + j*(Nx) + i) - num_of_evals_total

  IF (num_of_evals .GE. 0) THEN 
     DO count = 1,num_of_evals+1
        rbla = ran2(jdum)
     ENDDO
  ELSE
     jdum = idum_init
     DO count = 1, k*(Nx*Ny) + j*(Nx) + i +1
        rbla = ran2(jdum)
     ENDDO
  ENDIF
  num_of_evals_total = k*(Nx*Ny) + j*(Nx) + i + 1
  decomp_independent_random = rbla 

END FUNCTION decomp_independent_random
  

FUNCTION ran2(idum)
  INTEGER(kind=ki) :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  REAL(kind=kr)    :: ran2,AM,EPS,RNMX
  PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.D0/IM1,IMM1=IM1-1, &
     & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
     & NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  INTEGER(kind=ki) :: idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  IF (idum.LE.0) THEN
     idum=MAX(-idum,1)
     idum2=idum
     DO j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        IF (idum.LT.0) idum=idum+IM1
        IF (j.LE.NTAB) iv(j)=idum
     ENDDO
     iy=iv(1)
  ENDIF
  k=idum/IQ1
  idum=IA1*(idum-k*IQ1)-k*IR1
  IF (idum.LT.0) idum=idum+IM1
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2
  IF (idum2.LT.0) idum2=idum2+IM2
  j=1+iy/NDIV
  iy=iv(j)-idum2
  iv(j)=idum
  IF(iy.LT.1)iy=iy+IMM1
  ran2=MIN(AM*iy,RNMX)
  RETURN
END FUNCTION ran2
!  (C) Copr. 1986-92 Numerical Recipes Software 0!5,.
  
