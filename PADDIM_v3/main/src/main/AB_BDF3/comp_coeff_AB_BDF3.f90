! This subroutine computes the coeddicients in the AB/BDF3 time stepping scheme
SUBROUTINE comp_coeff_AB_BDF3(dt,dt1,dt2)
  USE defprecision_module
  IMPLICIT NONE
  REAL(kind=kr) :: dt,dt1,dt2
  REAL(kind=kr) :: delta1,delta2

  delta1 = dt1 / dt 
  delta2 = dt2 / dt

  ta0 = 1._kr + 1._kr / (1._kr + delta1) + 1._kr / (1._kr + delta1 + delta2) 
  ta0 = ta0 / dt
  ta1 = - ((1._kr + delta1) * (1._kr + delta1 + delta2)) / (delta1*(delta1 + delta2))
  ta1 = ta1 / dt
  ta2 =   (1._kr + delta1 + delta2) / (delta1 * delta2 * (1._kr + delta1))
  ta2 = ta2 / dt
  ta3 = - (1._kr + delta1) / ( delta2 * (delta1 + delta2) * (1._kr + delta1 + delta2) )
  ta3 = ta3 / dt 

  tb1 =   ((1._kr + delta1) * (1._kr + delta1 + delta2)) / (delta1*(delta1 + delta2))
  tb2 = - (1._kr + delta1 + delta2) / (delta1 * delta2)
  tb3 =   (1._kr + delta1) /  ( delta2 * (delta1 + delta2) )

END SUBROUTINE comp_coeff_AB_BDF3
