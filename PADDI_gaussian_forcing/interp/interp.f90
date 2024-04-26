!-----------------------------------!
!                                   !
!   File    :   interp.f90          !
!                                   !
!   Author  :   Dante Buhl          !
!                                   !
!   Date    :   Oct. 29, 2023       !
!                                   !
!-----------------------------------!

module interp

integer, parameter :: dp=selected_real_kind(15)
real(8), parameter :: pi = acos(-1.0_dp)
