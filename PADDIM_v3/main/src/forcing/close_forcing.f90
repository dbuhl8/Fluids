! File: close_forcing.f90
! Author: Dante Buhl
! Purpose: deallocate forcing arrays and write restart files for forcing

subroutine close_forcing

    use message_passing_module, ONLY: myid

    implicit none

    integer :: i

#ifdef STOCHASTIC_FORCING
    !write restart files
    if(forcingCPU) then
        open(69, file="forcing_data/restartforcing"//trim(str(myid))//".dat")
            do i = cTimeindex-window_pts+1, cTimeindex
                write(69, dataformat) gpTimevals(i), gpForcingVals(i, :)
            end do
        close(69)
        deallocate(gpForcingVals)
        deallocate(gpTimeVals)
        deallocate(interpSlopeVals)
    end if
    
#endif

    deallocate(waveNumMap)
    deallocate(force_real)
    deallocate(force_spec)
    

end subroutine close_forcing
