SUBROUTINE open_files
  USE defprecision_module
  USE parameter_module, ONLY : B_therm,B_comp,C_Lorentz,                            &
                             & D_visc,D_therm,D_comp,D_mag,                         &
                             & Nx,Ny,Nz,Lmax,Mmax,Nmax,                             &
                             & Gammax,Gammay,Gammaz,CFL_safety_fac,dt_max,          &
                             & dt_initial, S_therm, S_comp, R, Theta
  USE message_passing_module, ONLY : myid,nprocs1,nprocs2
  USE pnetCDF_IO_module, ONLY: pn_open_dump,pn_open_simdat_file,                    &
         & pn_open_xyslice_simdat_file,pn_open_xzslice_simdat_file,pn_open_yzslice_simdat_file
  IMPLICIT NONE
  INTEGER (kind=ki)            :: jcopen,jcwinfo2
  INTEGER (kind=ki),PARAMETER  :: ninfo=11
  REAL (kind=kr_jc)            :: finfo(ninfo)
  CHARACTER (LEN=20)           :: cinfo(ninfo)
  CHARACTER (LEN=100)          :: jcname
  CHARACTER (LEN=2)            :: acc
  INTEGER (kind=ki)            :: i,ierr
  INTEGER                      :: date_time(8)
  CHARACTER (LEN=12)           :: real_clock(3)
#ifdef TWO_DIMENSIONAL
  REAL(kind=kr_jc)             :: grid(2*nx*nz)
  INTEGER(kind=ki)             :: k,icwrite
#endif

! Only root process has to open files 
  IF (myid.EQ.0) THEN 
     !
     ! get date and time
     ! -----------------
     CALL DATE_AND_TIME(real_clock(1),real_clock(2),real_clock(3),date_time)
     !
     ! set output file units
     ! ---------------------
     DO i=1,nout
        uout(i)=80+i
     ENDDO
     !
     ! open output files 
     ! -----------------
     DO i=1,nout
        OPEN(UNIT=uout(i),FILE=outfile(i),ACTION="WRITE",POSITION="APPEND")
        WRITE(uout(i),'(a)')       "# =========================================================="
        WRITE(uout(i),'(a)')       "# ===            Double Diffusive Convection             ==="
        WRITE(uout(i),'(a)')       "# ==            Written by Stephan Stellmach              =="
        WRITE(uout(i),'(a)')       "# ==    Contact stellma@uni-muenster.de for questions     =="
        WRITE(uout(i),'(a)')       "# =========================================================="
        WRITE(uout(i),'(a)')       "#"
        WRITE(uout(i),'(a)')       "# Execution date and time:"
        WRITE(uout(i),'(a)')       "# ------------------------"
        WRITE(uout(i),'(a,i4,a,i2,a,i2,a,i2,a,i2,a,i2 )')                                  &
                                "# Run started: year:",date_time(1),                    &
&                               ", month:",date_time(2),", day:",date_time(3),          &
&                               ", Time:",date_time(5),":",date_time(6),":",date_time(7)
        WRITE(uout(i),'(a)')       "#"
        WRITE(uout(i),'(a)')       "# Parameter values:"
        WRITE(uout(i),'(a)')       "# -----------------"
        WRITE(uout(i),'(a,E20.7)') "# Thermal buoyancy coefficient        :",B_therm
        WRITE(uout(i),'(a,E20.7)') "# Compositional buoyancy coefficient  :",B_comp
        WRITE(uout(i),'(a,E20.7)') "# Lorentz force coefficient           :",C_Lorentz ! PH
        WRITE(uout(i),'(a,E20.7)') "# Viscous diffusion parameter         :",D_visc
        WRITE(uout(i),'(a,E20.7)') "# Thermal diffusion parameter         :",D_therm
        WRITE(uout(i),'(a,E20.7)') "# Compositional diffusion parameter   :",D_comp
        WRITE(uout(i),'(a,E20.7)') "# Magnetic diffusion parameter        :",D_mag ! PH
        WRITE(uout(i),'(a,E20.7)') "# Thermal background stratification parameter       :",S_therm
        WRITE(uout(i),'(a,E20.7)') "# Compositional background stratification parameter :",S_comp
        WRITE(uout(i),'(a,E20.7)') "# Rotational parameter :", R
        WRITE(uout(i),'(a,E20.7)') "# Angle Rotational Axis and Gravity parameter :",Theta 
        WRITE(uout(i),'(a,E20.7)') "# Length of box            :",Gammax
        WRITE(uout(i),'(a,E20.7)') "# Width of box             :",Gammay
        WRITE(uout(i),'(a,E20.7)') "# Height of box            :",Gammaz
        WRITE(uout(i),'(a,E20.7)') "# CFL_safety_factor 0",CFL_safety_fac
        WRITE(uout(i),'(a,E20.7)') "# Maximum allowed time step length =",dt_max
        WRITE(uout(i),'(a,E20.7)') "# Initial time step length =",dt_initial
        WRITE(uout(i),'(a,I6)')    "# Lmax = :",Lmax
        WRITE(uout(i),'(a,I6)')    "# Mmax = :",Mmax
        WRITE(uout(i),'(a,I6)')    "# Nmax = :",Nmax
        WRITE(uout(i),'(a,I6)')    "# Nx = :",Nx
        WRITE(uout(i),'(a,I6)')    "# Ny = :",Ny
        WRITE(uout(i),'(a,I6)')    "# Nz = :",Nz
        WRITE(uout(i),'(a,I6)')    "# number_of_tasks_1st_transpose = :",nprocs1
        WRITE(uout(i),'(a,I6)')    "# number_of_tasks_2nd_transpose = :",nprocs2
        WRITE(uout(i),'(a)')       "#"
        WRITE(uout(i),'(a)')       "# Compilation info - Compiled with following flags defined"
#ifdef TEMPERATURE_FIELD
        WRITE(uout(i),'(a)')       "#  -TEMPERATURE_FIELD"
#endif
#ifdef CHEMICAL_FIELD
        WRITE(uout(i),'(a)')       "#  -CHEMICAL_FIELD"
#endif
#ifdef MAGNETIC
        WRITE(uout(i),'(a)')       "#  -MAGNETIC"
#endif
#ifdef TWO_DIMENSIONAL
        WRITE(uout(i),'(a)')       "#  -TWO_DIMENSIONAL"
#endif
#ifdef AB_BDF3
        WRITE(uout(i),'(a)')       "#  -AB_BDF3"
#endif
#ifdef SINGLE_PRECISION
        WRITE(uout(i),'(a)')       "#  -SINGLE_PRECISION"
#endif
#ifdef DOUBLE_PRECISION
        WRITE(uout(i),'(a)')       "#  -DOUBLE_PRECISION"
#endif
        WRITE(uout(i),'(a)')       "#"
        WRITE(uout(i),'(a)')       "# =========================================================="
        WRITE(uout(i),*)
     ENDDO
     !
     ! open jc_file
     ! ------------
     IF (write_compressed_fields) THEN 
        acc = "w"
        jcname=jc_out
#ifdef TWO_DIMENSIONAL
        ierr = jcopen(1,JC2D,jcname,acc)
#else
        ierr = jcopen(1,JC3D,jcname,acc)
#endif
        cinfo(1)= "Thermal buoyancy coefficient"
        finfo(1)= B_therm
        cinfo(2)= "Compositional buoyancy coefficient"
        finfo(2)= B_comp
        cinfo(3)= "Lorentz force coefficient" ! PH
        finfo(3)= C_Lorentz ! PH
        cinfo(4)= "Viscous diffusion parameter"
        finfo(4)= D_visc
        cinfo(5)= "Thermal diffusion parameter"
        finfo(5)= D_therm
        cinfo(6)= "Compositional diffusion parameter"
        finfo(6)= D_comp
        cinfo(7)= "Magnetic diffusion parameter" ! PH
        finfo(7)= D_mag ! PH
        cinfo(8)= "Thermal background stratification parameter"
        finfo(8)= S_therm
        cinfo(9)= "Compositional background stratification parameter"
        finfo(9)= S_comp
        cinfo(10)= "Rotation parameter"
        finfo(10)= R
        cinfo(11)= "Angle of Rotational Axis and  Gravity parameter"
        finfo(11)= Theta
       
#ifdef TWO_DIMENSIONAL
        number_of_jc_fields = 3
#else
        number_of_jc_fields = 6
#endif
#ifdef TEMPERATURE_FIELD
        number_of_jc_fields = number_of_jc_fields + 1 
#endif
#ifdef CHEMICAL_FIELD
        number_of_jc_fields = number_of_jc_fields + 1 
#endif
#ifdef MAGNETIC
        number_of_jc_fields = number_of_jc_fields + 3
#endif
#ifdef TWO_DIMENSIONAL
        ierr = jcwinfo2(1,nx,nz,1,number_of_jc_fields,ninfo,finfo,cinfo)
        ! Write grid in 2D
        ! CAUTION: Since the tool "rasterfari" used for visualization crashes
        !          for box heights different from 1, we have to take that into
        !          account.This is not nice of course,so beware!
        DO k=0,Nz-1
           DO i=0,Nx-1
              ! x-coordinate
              grid(i*nz+(k+1)) = (REAL(i,kr_jc) / REAL(nx-1,kr_jc)) * gammax/gammaz
              ! z-coordinate 
              grid(i*nz + (k+1) + nx*nz) = (REAL(k,kr_jc) / REAL(nz-1,kr_jc))
           END DO
        END DO
        WRITE(*,*) "1"
        ierr = icwrite(1,0,0._kr_jc,0._kr_jc,grid,2*nx*nz)
        WRITE(*,*) "2"
#else
        ierr = jcwinfo2(1,nx,ny,nz,number_of_jc_fields,ninfo,finfo,cinfo)
#endif
     ENDIF
        
  ENDIF

  ! Open netCDF Dump File
  CALL pn_open_dump
  ! Open netCDF simdat file
  CALL pn_open_simdat_file
#ifndef TWO_DIMENSIONAL
  CALL pn_open_xyslice_simdat_file
  CALL pn_open_xzslice_simdat_file
  CALL pn_open_yzslice_simdat_file
#endif

  ! print*,"Hallo"  

END SUBROUTINE open_files
