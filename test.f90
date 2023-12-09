program test
integer :: i
do i=1, 100
    open(11, file='Output'//trim(str(i))//'.txt')
    write (11, *) i
    close (11)
end do

contains 
character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

end program test



NAMELIST /input_values/ Restart_from_dumped_data
NAMELIST /input_values/ Restart_from_dumped_data_hydro
NAMELIST /input_values/ Restart_from_netCDF_output_file
NAMELIST /input_values/ max_degree_of_x_fourier_modes
NAMELIST /input_values/ max_degree_of_y_fourier_modes
NAMELIST /input_values/ max_degree_of_z_fourier_modes
NAMELIST /input_values/ dealiasing
NAMELIST /input_values/ Thermal_buoyancy_param
NAMELIST /input_values/ Compositional_buoyancy_param
NAMELIST /input_values/ Lorentz_force_coeff !PH
NAMELIST /input_values/ Viscous_diffusion_coeff
NAMELIST /input_values/ Thermal_diffusion_coeff
NAMELIST /input_values/ Compositional_diffusion_coeff
NAMELIST /input_values/ Magnetic_diffusion_coeff !PH
NAMELIST /input_values/ Thermal_stratif_param
NAMELIST /input_values/ Compositional_stratif_param
NAMELIST /input_values/ Rotational_param
NAMELIST /input_values/ Angle_rot_axis_gravity
NAMELIST /input_values/ x_extent_of_the_box
NAMELIST /input_values/ y_extent_of_the_box
NAMELIST /input_values/ z_extent_of_the_box
NAMELIST /input_values/ initial_time_step_length
NAMELIST /input_values/ maximum_time_step_length
NAMELIST /input_values/ CFL_safety_factor
NAMELIST /input_values/ number_of_time_steps
NAMELIST /input_values/ save_state_every_nth_timestep
NAMELIST /input_values/ save_state_netcdf_ev_nth_step
NAMELIST /input_values/ save_slice_netcdf_ev_nth_step
NAMELIST /input_values/ comp_diagno_every_nth_timestep
NAMELIST /input_values/ restart_info_every_nth_timestep
NAMELIST /input_values/ write_spec_every_nth_timestep


