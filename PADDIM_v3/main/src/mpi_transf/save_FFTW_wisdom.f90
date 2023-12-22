SUBROUTINE save_FFTW_wisdom
!+--------------------------------------------------------------------+
!| This subroutine uses the wisdom mechanism of FFTW to write         |
!| information collected during the FFTW-plannung phase to disk.      |
!+--------------------------------------------------------------------+
  IMPLICIT NONE
  IF (use_FFTW_wisdom_file) THEN 
     OPEN(UNIT=uwisdom,FILE=FFTW_wisdom_file)
     REWIND(UNIT=uwisdom)
     CALL export_wisdom_to_file(uwisdom)
     CLOSE(uwisdom)
  ENDIF
END SUBROUTINE save_FFTW_wisdom
