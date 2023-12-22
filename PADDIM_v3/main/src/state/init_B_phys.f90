SUBROUTINE init_B_phys(B)
  USE defprecision_module
  USE parameter_module, ONLY: Nx,Ny,Nz,Gammax,Gammay,Gammaz,pi,alpha,beta,gamma
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys
  USE message_passing_module, ONLY: myid
  IMPLICIT NONE
  TYPE(field) :: B
  INTEGER(kind=ki) :: i,j,k
  REAL(kind=kr) :: xc,yc,zc,dx,dy,dz, xmid, ymid, zmid, Bscale

  dx = Gammax / REAL(Nx, kind=kr)
  dy = Gammay / REAL(Ny, kind=kr)
  dz = Gammaz / REAL(Nz, kind=kr)
  
  Bscale = 1._kr
  xmid = 0.5_kr*Gammax
  ymid = 0.5_kr*Gammay
  zmid = 0.5_kr*Gammaz
  
  IF (myid .EQ. 0) THEN
      PRINT *, "Gammaz=", Gammaz, "Nz=", Nz, "dz=", Gammaz/REAL(Nz,kind=kr), "zmid=", zmid
  ENDIF
      
  DO k=mysz_phys,myez_phys
     zc = k*dz
     DO j=mysy_phys,myey_phys
        yc = j*dy
        DO i=0,Nx-1
           xc = i*dx
#ifdef TWO_DIMENSIONAL
           ! Initial B field for 2D case 
           B%phys(i,j,k,vec_x) =  0._kr
           B%phys(i,j,k,vec_z) =  0._kr
#else
           ! Initial B field for 3D case 
           B%phys(i,j,k,vec_x) =  1._kr !1._kr/SQRT(2._kr)
           B%phys(i,j,k,vec_y) =  0._kr
           B%phys(i,j,k,vec_z) =  0._kr !1._kr/SQRT(2._kr)
#endif
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE init_B_phys
