! The follwing subroutine transposes a three dimensional 
! array in(:,:,:) which is decomposed along the x1-direction to an 
! array out(:,:,:) which is decomposed along the x2-direction. 
! The order in which the different directions are stored changes 
! from x2:x3:x1 -> x1:x2:x3.  
SUBROUTINE transp1d_231_1dec_to_123_2dec(slice2,slice1,info)
  USE defprecision_module
  IMPLICIT NONE
  TYPE(transpose_info_1d_decomp) :: info
  COMPLEX(kind=kr)               :: slice1(:,info%mys2:,:)
  COMPLEX(kind=kr)               :: slice2(:,:,info%mys1:)
  COMPLEX(kind=kr),ALLOCATABLE   :: work(:,:,:)
  INTEGER(kind=ki)               :: k

! allocate temporary memory 
  ALLOCATE(work(info%mys1:info%mye1,info%M3,info%N2))

! Transpose from x2:x3:x1 -> x1:x3:x2 in local memory
  DO k=1,info%M3
     work(:,k,:) = TRANSPOSE(slice2(:,k,:))
  ENDDO
  ! second possible way using RESHAPE: seems to be slower with Intel Fortran 
  ! work = reshape(slice2,(/ info%mye1-info%mys1+1,info%M3,infoa%N2/), &
  ! &               order = (/3,2,1/)                                  )
 
! Transpose from x1:x3:x2 with x1-decomp. -> x1:x2:x3 with x2-decomp. 
  call transp1d_132_1dec_to_123_2dec(work,slice1,info)

! deallocate temporary memory
  deallocate(work)

END SUBROUTINE transp1d_231_1dec_to_123_2dec

