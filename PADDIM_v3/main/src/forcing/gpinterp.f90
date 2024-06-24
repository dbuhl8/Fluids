! File: interp.f90
! Author: Dante Buhl
! Description: computes a linear regression of the forcing columns in GP

function gpinterp(t, i, j) result(output)

    use defprecision_module   
 
    implicit none

    integer(kind=ki) :: i, j, k, cfi
    real(kind=kr) :: t
    real(kind=kr), dimension(2) :: output

    !This updates currentTimeIndex
    if (t .ge. gpTimeVals(cTimeIndex)) then
      do k = cTimeIndex+1, numGProws 
        if (t .lt. gpTimeVals(k)) then 
          cTimeIndex= k
          exit
        end if
      end do
    end if


    !build the linear interpolation
    cfi = 2*waveNumMap(i, j)
   
    output(1) = ((gpForcingVals(cTimeIndex,cfi-1)-gpForcingVals(cTimeIndex-1,cfi-1))/&
      tstep) * (t - gpTimeVals(cTimeIndex-1))&
      + gpForcingVals(cTimeIndex-1, cfi-1)
    output(2) = ((gpForcingVals(cTimeIndex,cfi)-gpForcingVals(cTimeIndex-1,cfi))/&
      tstep) * (t - gpTimeVals(cTimeIndex-1))&
      + gpForcingVals(cTimeIndex-1, cfi)

    ! old method using a slope array (got rid of this to improve memory usage)
    !output(1)  = interpSlopeVals(cTimeIndex-1, cfi-1) * (t - gpTimeVals(cTimeIndex-1)) + &
      !gpForcingVals(cTimeIndex-1, cfi-1)
    !output(2)  = interpSlopeVals(cTimeIndex-1, cfi) * (t - gpTimeVals(cTimeIndex-1)) + &
      !gpForcingVals(cTimeIndex-1, cfi)

end function gpinterp
