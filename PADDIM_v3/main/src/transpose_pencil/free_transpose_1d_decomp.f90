SUBROUTINE free_transpose_1d_decomp(info)
  IMPLICIT NONE
  TYPE(transpose_info_1d_decomp) :: info
  
  ! deallocate memory 

  DEALLOCATE(info%s1)
  DEALLOCATE(info%e1)
  DEALLOCATE(info%s2)
  DEALLOCATE(info%e2)

  DEALLOCATE(info%sendcounts)
  DEALLOCATE(info%sdispl)
  DEALLOCATE(info%recvcounts)
  DEALLOCATE(info%rdispl)

END SUBROUTINE free_transpose_1d_decomp
