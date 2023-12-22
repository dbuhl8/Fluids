! The code was written to run efficiently in 2D and 3D. 
! Both cases are largely equivalent. However, in 3D, both a vectorfield and 
! its curl have 3 components, whereas in 2D, the vectorfield has 2 and the 
! curl only 1 nonzero components. We therefore define the following quanteties 
! to make indexing throughout the code as easy as possible:
module defs_2D_3D_module
  use defprecision_module
  implicit none 

#ifdef TWO_DIMENSIONAL
  INTEGER(kind=ki), PARAMETER :: dim_vec = 2   ! dimensions of vector field
  INTEGER(kind=ki), PARAMETER :: dim_curl = 1  ! dimensions of the curl 
  ! velocity
  INTEGER(kind=ki), PARAMETER :: vec_x = 1 ! x-component of a vector field
  INTEGER(kind=ki), PARAMETER :: vec_z = 2 ! z-component of a vector field
  ! vorticity
  INTEGER(kind=ki), PARAMETER :: curl_y = 1 !y-component of the curl
#else 
  INTEGER(kind=ki), PARAMETER :: dim_vec = 3
  INTEGER(kind=ki), PARAMETER :: dim_curl = 3
  ! velocity
  INTEGER(kind=ki), PARAMETER :: vec_x = 1 
  INTEGER(kind=ki), PARAMETER :: vec_y = 2 
  INTEGER(kind=ki), PARAMETER :: vec_z = 3 
  ! vorticity
  INTEGER(kind=ki), PARAMETER :: curl_x = 1 
  INTEGER(kind=ki), PARAMETER :: curl_y = 2 
  INTEGER(kind=ki), PARAMETER :: curl_z = 3 
#endif
end module defs_2D_3D_module
