!+--------------------------------------------------------------------+
!| The following function computes the minumum and maximum of         |
!|         T,C,u_x,u_y,u_z,VORT_x,VORT_y,VORT_z, Bx, By, Bz           |
!| as well as the maxima of                                           |
!|                | u | , | VORT | , and | B |.                       |
!| Caution: It is assumed that this routine is not called at every    |
!|          time step so that it does not need to be optimized.       |
!|          In fact it is likely to be highly inefficient...          |
!+--------------------------------------------------------------------+
SUBROUTINE compute_uTCB_minmax(Temp_min,Temp_max,Chem_min,Chem_max,              &
           &                  u_min,u_max,VORT_min,VORT_max,B_min,B_max,         &
           &                  u_max_abs,VORT_max_abs,B_max_abs,                  &
           &                  u,Temp,Chem,B) ! PH
  USE defprecision_module
  USE state_module,ONLY : velocity,buoyancy,field ! PH
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  TYPE(velocity) :: u
  TYPE(field)    :: B ! PH
  TYPE(buoyancy) :: Temp,Chem
  REAL(kind=kr) :: Temp_min,Temp_max,Chem_min,Chem_max,           &
       &           u_min(3),u_max(3),VORT_min(3),VORT_max(3),     &
       &           B_min(3), B_max(3), u_max_abs,VORT_max_abs, B_max_abs ! PH
  REAL(kind=kr) :: min_local,max_local,nmin_local(3),nmax_local(3)
  INTEGER(kind=ki) :: ierr

  ! Temperature fluctuation
#ifdef TEMPERATURE_FIELD
  min_local = MINVAL(Temp%phys)
  max_local = MAXVAL(Temp%phys) 
  CALL MPI_REDUCE(min_local,Temp_min,1,PM_MPI_FLOAT_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(max_local,Temp_max,1,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)
#else
  Temp_min = 0._kr
  Temp_max = 0._kr
#endif

  ! Compositional fluctuation

#ifdef CHEMICAL_FIELD
  min_local = MINVAL(Chem%phys)
  max_local = MAXVAL(Chem%phys) 
  CALL MPI_REDUCE(min_local,Chem_min,1,PM_MPI_FLOAT_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(max_local,Chem_max,1,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)
#else 
  Chem_min = 0._kr
  Chem_max = 0._kr
#endif

  ! velocity components
  nmin_local(1) = MINVAL(u%phys(:,:,:,vec_x))
  nmax_local(1) = MAXVAL(u%phys(:,:,:,vec_x))
#ifndef TWO_DIMENSIONAL
  nmin_local(2) = MINVAL(u%phys(:,:,:,vec_y))
  nmax_local(2) = MAXVAL(u%phys(:,:,:,vec_y))
#else
  nmin_local(2) = 0._kr 
  nmax_local(2) = 0._kr
#endif
  nmin_local(3) = MINVAL(u%phys(:,:,:,vec_z))
  nmax_local(3) = MAXVAL(u%phys(:,:,:,vec_z))

  CALL MPI_REDUCE(nmin_local,u_min,3,PM_MPI_FLOAT_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(nmax_local,u_max,3,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)

  ! vorticity components
  nmin_local(2) = MINVAL(u%curl(:,:,:,curl_y))
  nmax_local(2) = MAXVAL(u%curl(:,:,:,curl_y))
#ifndef TWO_DIMENSIONAL
  nmin_local(1) = MINVAL(u%curl(:,:,:,curl_x))
  nmax_local(1) = MAXVAL(u%curl(:,:,:,curl_x))
  nmin_local(3) = MINVAL(u%curl(:,:,:,curl_z))
  nmax_local(3) = MAXVAL(u%curl(:,:,:,curl_z))
#else 
  nmin_local(1) = 0._kr
  nmax_local(1) = 0._kr
  nmin_local(3) = 0._kr
  nmax_local(3) = 0._kr
#endif
  CALL MPI_REDUCE(nmin_local,VORT_min,3,PM_MPI_FLOAT_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(nmax_local,VORT_max,3,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)
 
#ifdef MAGNETIC  
  ! B field components
  nmin_local(1) = MINVAL(B%phys(:,:,:,vec_x))
  nmax_local(1) = MAXVAL(B%phys(:,:,:,vec_x))
#ifndef TWO_DIMENSIONAL
  nmin_local(2) = MINVAL(B%phys(:,:,:,vec_y))
  nmax_local(2) = MAXVAL(B%phys(:,:,:,vec_y))
#else
  nmin_local(2) = 0._kr 
  nmax_local(2) = 0._kr
#endif
  nmin_local(3) = MINVAL(B%phys(:,:,:,vec_z))
  nmax_local(3) = MAXVAL(B%phys(:,:,:,vec_z))  

  CALL MPI_REDUCE(nmin_local,B_min,3,PM_MPI_FLOAT_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(nmax_local,B_max,3,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)
#else
  B_min = 0._kr
  B_max = 0._kr
#endif


  ! maximum absolute value of velocity
#ifdef TWO_DIMENSIONAL
  max_local = MAXVAL( SQRT( u%phys(:,:,:,vec_x)**2 + u%phys(:,:,:,vec_z)**2 ) )
#else
  max_local = MAXVAL( SQRT( u%phys(:,:,:,vec_x)**2 &
              &           + u%phys(:,:,:,vec_y)**2 &
              &           + u%phys(:,:,:,vec_z)**2 ) )
#endif
  CALL MPI_REDUCE(max_local,u_max_abs,1,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)

  ! maximum absolute value of vorticity
#ifdef TWO_DIMENSIONAL
  max_local = MAXVAL( SQRT( u%curl(:,:,:,curl_y)**2 ) )
#else
  max_local = MAXVAL( SQRT( u%curl(:,:,:,curl_x)**2 &
              &           + u%curl(:,:,:,curl_y)**2 &
              &           + u%curl(:,:,:,curl_z)**2 ) )
#endif
  CALL MPI_REDUCE(max_local,VORT_max_abs,1,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  
#ifdef MAGNETIC
  ! maximum absolute value of B field
#ifdef TWO_DIMENSIONAL
  max_local = MAXVAL( SQRT( B%phys(:,:,:,vec_x)**2 + B%phys(:,:,:,vec_z)**2 ) )
#else
  max_local = MAXVAL( SQRT( B%phys(:,:,:,vec_x)**2 &
              &           + B%phys(:,:,:,vec_y)**2 &
              &           + B%phys(:,:,:,vec_z)**2 ) )
#endif
  CALL MPI_REDUCE(max_local,B_max_abs,1,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)
#else
  B_max_abs = 0._kr
#endif

END SUBROUTINE compute_uTCB_minmax
