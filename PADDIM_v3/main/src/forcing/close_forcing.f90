! File: close_forcing.f90
! Author: Dante Buhl
! Purpose: deallocate forcing arrays

subroutine close_forcing
  implicit none

  deallocate(waveNumMap)
  deallocate(force_real)
  deallocate(force_spec)
end subroutine close_forcing
