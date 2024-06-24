! File: close_forcing.f90
! Author: Dante Buhl
! Purpose: deallocate forcing arrays and write restart files for forcing

subroutine write_forcing_dump
  use message_passing_module, ONLY: myid
  implicit none
  integer :: i

  !write restart files
  if(forcingCPU) then
    open(100, file="forcing_data/"//fdump_out_file//"_forcing"//&
      trim(str(myid))//".dat")
      ! this needs to be changed in order to write out entire forcing array
      write(100, headerFormatInt) "Number of GP Columns :", numGPcolumns
      write(100, headerFormatInt) "Number of GP Rows    :", numGProws
      write(100, headerFormatInt) "Random Num Index     :", fidum
      do i = 1, numGProws
        write(100, dataformat) gpTimevals(i), gpForcingVals(i, :)
      end do
    close(100)
  end if
end subroutine write_forcing_dump
