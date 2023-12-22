SUBROUTINE adapt_dt_AB_BDF3(u,dt,dt1)
  USE defprecision_module
  USE state_module, ONLY: velocity
  USE parameter_module, ONLY: dt_max,CFL_safety_fac
  IMPLICIT NONE
  TYPE(velocity) :: u
  REAL(kind=kr)  :: dt,dt1
  REAL(kind=kr)  :: dt_cfl
  
  dt_cfl = CFL_safety_fac*CFL_AB_BDF3(u)

  IF (dt_cfl .GT. 1.05*dt1) THEN !allow to increase the time step by 5% maximum
     dt = MIN(1.05*dt1,dt_max)
  ELSE
     dt = MIN(dt_cfl,dt_max)
  ENDIF

END SUBROUTINE adapt_dt_AB_BDF3

