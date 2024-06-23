!Routine: compute_weighted_u_z
!Author: Dante Buhl
!Purpose: To compute rms values for u_z weighted with the vertical enstrophy as an indicator for turbulence.
!Conceptualized by Pascale Garaud


subroutine compute_weighted_u_z(u, VORTZ_rms, u_turb, u_lam)
  ! args needed u, Vortz_rms (already computed by compute_uTCB_rms), u_turb, u_lam)
                
  ! This is an excellent indicator for the local small scale turbulence. 
 
  use defprecision_module
  use defs_2D_3D_module, ONLY: vec_z, curl_z
  use state_module, ONLY: velocity
  use parameter_module, ONLY: Nx
  use mpi_transf_module, ONLY:  mysy_phys,myey_phys,mysz_phys,myez_phys
  use MPI

  implicit none

  real(kind=kr)  :: u_turb, u_lam
  real(kind=kr), intent(in) :: VORTZ_rms
  type(velocity) :: u
  real(kind=kr), allocatable :: work(:, :, :)
  real(kind=kr) :: u_z_dot_w_z_rms, u_z_over_w_z_rms, one_over_w_z_rms

  allocate(work(0:Nx-1, mysy_phys:myey_phys, mysz_phys:myez_phys))

    ! compute u_z * w_z rms
    work = u%phys(:, :, :, vec_z) * u%curl(:, :, :, curl_z)
    u_z_dot_w_z_rms = rms(work)

    ! compute 1 / w_z rms
    work = 1. / u%curl(:, :, :, curl_z)

    !work2 = 1. - u%curl(:, :, :, curl_z)
    !call geosum(work2, 15)

    ! this statement reduces numerical precision errors caused by division by 0 or near zero.
    !where (abs(u%curl(:, :, :, curl_z)) .lt. 1.0_kr) work = work2

    one_over_w_z_rms = rms(work)

    ! compute u_z / w_z rms
    work = work * u%phys(:, :, :, vec_z)
    u_z_over_w_z_rms = rms(work)

    
    ! calculate u_turb, u_lam
    u_turb = u_z_dot_w_z_rms/VORTZ_rms
    u_lam  = u_z_over_w_z_rms/one_over_w_z_rms

    deallocate(work)

end subroutine compute_weighted_u_z
