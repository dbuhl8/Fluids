! The follwing subroutine transposes a three dimensional 
! array in(:,:,:) which is decomposed along the x2-direction to an 
! array out(:,:,:) which is decomposed along the x1-direction. 
! The order in which the different directions are stored changes 
! from x1:x2:x3 -> x2:x3:x1.  
SUBROUTINE transp1d_123_2dec_to_231_1dec(slice1,slice2,info)
  USE defprecision_module
  IMPLICIT NONE
  TYPE(transpose_info_1d_decomp) :: info
  COMPLEX(kind=kr)               :: slice1(:,info%mys2:,:)
  COMPLEX(kind=kr)               :: slice2(:,:,info%mys1:)
  COMPLEX(kind=kr),ALLOCATABLE   :: work(:,:,:)
  INTEGER(kind=ki)               :: k

! allocate temporary memory 
  ALLOCATE(work(info%mys1:info%mye1,info%M3,info%N2))

! Transpose from x1:x2:x3 with x2-decomposition -> x1:x3:x2 with x1-decomp.
  call transp1d_123_2dec_to_132_1dec(slice1,work,info)

! finally, Transpose from x1:x3:x2 -> x2:x3:x1 in local memory
  DO k=1,info%M3   
     slice2(:,k,:) = TRANSPOSE(work(:,k,:))
  ENDDO
  ! RESHAPE may also be used: seems to be slower with Intel Fortran 
  ! slice2 = reshape(work,(/info%N2,info%M3,info%mye1-info%mys1+1/), &
  ! &                order=(/3,2,1/)                                  )

! deallocate temporary memory
  DEALLOCATE(work)

END SUBROUTINE transp1d_123_2dec_to_231_1dec
