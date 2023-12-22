SUBROUTINE get_FFTW_wisdom
!+--------------------------------------------------------------------+
!| This subroutine uses the wisdom mechanism of FFTW to read          |
!| information which can be used to speed up the FFTW plan computaion |
!| from disk.                                                         |
!+--------------------------------------------------------------------+
  IMPLICIT NONE
  INTEGER :: isuccess
  IF (use_FFTW_wisdom_file) THEN 
     OPEN(UNIT=uwisdom,FILE=FFTW_wisdom_file)
     REWIND(UNIT=uwisdom)
     CALL import_wisdom_from_file(isuccess,uwisdom)
     CLOSE(uwisdom)
     IF (isuccess.NE.1) THEN 
        WRITE(*,'(a,a,a)') "FFTW wisdom import from file ",       &
&                           TRIM(FFTW_wisdom_file)," not succussful"
     ELSE
        WRITE(*,'(a,a,a)') "FFTW wisdom import from file ",       &
&                           TRIM(FFTW_wisdom_file)," succussful"      
     ENDIF
  ENDIF
END SUBROUTINE get_FFTW_wisdom

