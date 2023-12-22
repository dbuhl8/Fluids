! File: close_stochastic_forcing.f90
! Author: Dante Buhl
! Purpose: deallocate forcing arrays

subroutine close_stochastic_forcing
    implicit none

    
    deallocate(waveNumMap)
    deallocate(gpForcingVals)
    deallocate(gpTimeVals)
    deallocate(interpSlopeVals)
    

end subroutine close_stochastic_forcing
