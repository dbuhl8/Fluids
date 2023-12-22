SUBROUTINE shift_time_pointers
  IMPLICIT NONE
  INTEGER(kind=ki) :: ihelp

  ihelp = ltime2
  ltime2 = ltime1
  ltime1 = ltime0 
  ltime0 = ihelp

  ihelp = rtime2
  rtime2 = rtime1
  rtime1 = ihelp

END SUBROUTINE shift_time_pointers
