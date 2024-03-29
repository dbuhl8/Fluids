! File: interp.f90
! Author: Dante Buhl
! Description: computes a linear regression of the forcing columns in GP

function gpinterp(t, i, j) result(output)
    
    implicit none

    integer(kind=ki) :: i, j, k, cForcingIndex
    real(kind=kr) :: t
    real(kind=kr), dimension(2) :: output


    !This updates currentTimeIndex
    if (t .ge. gpTimeVals(cTimeIndex)) then
        do k = cTimeIndex, numGProws 
            if (t .lt. gpTimeVals(k)) then 
                cTimeIndex= k
                exit
            end if
        end do
    end if


    !build the linear interpolation
    cForcingIndex = 2*waveNumMap(i, j)
   
    output(1)  = interpSlopeVals(cTimeIndex-1, cForcingIndex-1) * (t - gpTimeVals(cTimeIndex-1)) + &
                              gpForcingVals(cTimeIndex-1, cForcingIndex-1)
    output(2)  = interpSlopeVals(cTimeIndex-1, cForcingIndex) * (t - gpTimeVals(cTimeIndex-1)) + &
                              gpForcingVals(cTimeIndex-1, cForcingIndex)

end function gpinterp
