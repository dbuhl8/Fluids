SUBROUTINE decomp_1d( n, numtasks, myid, s, e )
  USE defprecision_module
  IMPLICIT NONE
  INTEGER (kind=ki) :: n, numtasks, myid, s, e
  INTEGER (kind=ki) :: nlocal
  INTEGER (kind=ki) :: deficit

  nlocal  = n / numtasks
  s       = myid * nlocal + 1
  deficit = MOD(n,numtasks)
  s       = s + MIN(myid,deficit)
  IF (myid .LT. deficit) THEN
     nlocal = nlocal + 1
  ENDIF
  e = s + nlocal - 1
  IF (e .GT. n .OR. myid .EQ. numtasks-1) e = n

END SUBROUTINE decomp_1d
