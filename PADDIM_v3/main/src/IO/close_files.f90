SUBROUTINE close_files
  USE defprecision_module
  USE message_passing_module, ONLY: myid
  USE pnetCDF_IO_module, ONLY: pn_close_dump,pn_close_simdat_file,pn_close_slice_simdat_file
  IMPLICIT NONE
  INTEGER (kind=ki) :: jcclose
  INTEGER (kind=ki) :: i,ierr

  ! only root process has to close files
  IF (myid.EQ.0) THEN
     ! close output files
     DO i=1,nout
        CLOSE(UNIT=uout(i))
     ENDDO
     ! close JC-Files
     IF (write_compressed_fields) THEN
        ierr=jcclose(1)
     ENDIF
  ENDIF

  ! close netCDF restart files
  CALL pn_close_dump
  CALL pn_close_simdat_file
#ifndef TWO_DIMENSIONAL
  CALL pn_close_slice_simdat_file
#endif

END SUBROUTINE close_files
