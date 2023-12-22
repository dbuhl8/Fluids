SUBROUTINE free_allocated_memory(u,Temp,Chem,B)
  USE defprecision_module
  USE state_module, ONLY: velocity,buoyancy,field,deallocate_uTCB
  USE mpi_transf_module, ONLY: free_transforms
  USE parameter_module, ONLY: deallocate_fft_storage_scheme
  IMPLICIT NONE
  TYPE(velocity) :: u
  TYPE(field)    :: B
  TYPE(buoyancy) :: Temp,Chem

  CALL free_transforms
  CALL deallocate_fft_storage_scheme
  CALL deallocate_uTCB(u,Temp,Chem,B)

END SUBROUTINE free_allocated_memory
