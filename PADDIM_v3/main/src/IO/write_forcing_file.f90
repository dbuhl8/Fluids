! File: close_forcing.f90
! Author: Dante Buhl
! Purpose: deallocate forcing arrays and write restart files for forcing

subroutine write_forcing

    use message_passing_module, ONLY: myid
    use forcing_module

    implicit none

    integer :: i

#ifdef STOCHASTIC_FORCING
    !write restart files
    if(forcingCPU) then
        open(69, file="forcing_data/restartforcing"//trim(str(myid))//".dat")
            ! need to adjust this to be the length of a dump span
            ! save (n_wrt_dump*dt)(gst) + window_pts 
            ! what might be done is to print the number of points to be read 
            ! into the next file at the very top. 
            do i = cTimeindex-window_pts+1, cTimeindex
                write(69, dataformat) gpTimevals(i), gpForcingVals(i, :)
            end do
        close(69)
    end if
    
#endif
    

end subroutine close_forcing
