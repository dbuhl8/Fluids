! File: interp.f90
! Author: Dante Buhl

function gpinterp(t, i, j) result(output)
    
    use message_passing_module, ONLY: myid

    implicit none

    integer(kind=ki) :: i, j, k, cForcingIndex
    real(kind=kr) :: t, output

    print *, "cpu "//trim(str(myid))//": entered gpinterp call"

    !This updates currentTimeIndex
    if (t .ge. gpTimeVals(cTimeIndex)) then
        do k = cTimeIndex, numGProws 
            if (t .lt. gpTimeVals(k)) then 
                cTimeIndex= k
                exit
            end if
        end do
    end if

    print *, "cpu "//trim(str(myid))//": updated time index correctly"

    !build the linear interpolation
    
    print *, "cpu "//trim(str(myid))//": found GP index for mode", i, j
    
    cForcingIndex = waveNumMap(i, j)
   
    print *, "cpu "//trim(str(myid))//": found GP index: ", cForcingIndex
     
    output  = interpSlopeVals(cTimeIndex, cForcingIndex) * (t - gpTimeVals(cTimeIndex-1)) + &
                              gpForcingVals(cTimeIndex-1, cForcingIndex)

    print *, "cpu "//trim(str(myid))//": computed interpolation"
  
    contains

    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

end function gpinterp
