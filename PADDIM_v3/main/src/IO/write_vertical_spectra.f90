SUBROUTINE write_vertical_spectra(u,Temp,Chem,B,istep,t)
  USE defprecision_module
  USE state_module, ONLY: buoyancy, velocity,field,ltime0 ! PH
  USE message_passing_module, ONLY : myid
  USE parameter_module, ONLY: kz,Nmax
  USE mpi_transf_module,  ONLY:  mysx_spec,myex_spec,mysy_spec,myey_spec
  IMPLICIT NONE
  TYPE(velocity)   :: u
  TYPE(field)      :: B ! PH
  TYPE(buoyancy)   :: Temp,Chem
  INTEGER(kind=ki) :: istep
  REAL(kind=kr)    :: t
  REAL(kind=kr), POINTER   :: Ener_spec_u(:,:), Ener_spec_B(:,:) ! PH
  REAL(Kind=kr), POINTER   :: Ener_spec_Temp(:),Ener_spec_Chem(:)
  REAL(kind=kr)    :: Ekin,ETemp,EChem,Emag ! PH
  INTEGER(kind=ki) :: n

  ! allocate arrays for the vertical energy spectra
  ALLOCATE(Ener_spec_u(0:2*Nmax-1,3))
  ALLOCATE(Ener_spec_B(0:2*Nmax-1,3)) ! PH
  ALLOCATE(Ener_spec_Temp(0:2*Nmax-1))
  ALLOCATE(Ener_spec_Chem(0:2*Nmax-1))

  ! compute spectra

  Ener_spec_u(:,1) = vertical_power_spectrum(u%spec(:,:,:,vec_x,ltime0))
#ifdef TWO_DIMENSIONAL
  Ener_spec_u(:,2) = 0._kr
#else
  Ener_spec_u(:,2) = vertical_power_spectrum(u%spec(:,:,:,vec_y,ltime0))
#endif
  Ener_spec_u(:,3) = vertical_power_spectrum(u%spec(:,:,:,vec_z,ltime0))
#ifdef TEMPERATURE_FIELD
  Ener_spec_Temp   = vertical_power_spectrum(Temp%spec(:,:,:,ltime0))
#else
  Ener_spec_Temp   = 0._kr
#endif
#ifdef CHEMICAL_FIELD
  Ener_spec_Chem   = vertical_power_spectrum(Chem%spec(:,:,:,ltime0))
#else
  Ener_spec_Chem   = 0._kr
#endif

#ifdef MAGNETIC
  Ener_spec_B(:,1) = vertical_power_spectrum(B%spec(:,:,:,vec_x,ltime0))
#ifdef TWO_DIMENSIONAL
  Ener_spec_B(:,2) = 0._kr
#else
  Ener_spec_B(:,2) = vertical_power_spectrum(B%spec(:,:,:,vec_y,ltime0))
#endif
  Ener_spec_B(:,3) = vertical_power_spectrum(B%spec(:,:,:,vec_z,ltime0))
#else
  Ener_spec_B   = 0._kr
#endif


  ! write spectrum
  ! PH{
  IF (myid.EQ.0) THEN 
     WRITE (uout(3),'(a,I8,a,E20.7)') '# Timstep =',istep,' time =',t
     Ekin=0._kr ! used for testing: see below
     Etemp=0._kr
     EChem=0._kr
     DO n=0,2*Nmax-1
        WRITE(uout(3),'(11E20.7)') kz(n),                      &
             & Ener_spec_u(n,1) , Ener_spec_u(n,2) , Ener_spec_u(n,3), &
             & Ener_spec_u(n,1) + Ener_spec_u(n,2) + Ener_spec_u(n,3), &
             & Ener_spec_B(n,1) , Ener_spec_B(n,2) , Ener_spec_B(n,3), &
             & Ener_spec_B(n,1) + Ener_spec_B(n,2) + Ener_spec_B(n,3), &
             & Ener_spec_Temp(n),Ener_spec_Chem(n)
        Ekin=Ekin+(Ener_spec_u(n,1) + Ener_spec_u(n,2) + Ener_spec_u(n,3))
        Emag=Emag+(Ener_spec_B(n,1) + Ener_spec_B(n,2) + Ener_spec_B(n,3))
        Etemp=Etemp+Ener_spec_Temp(n)
        Echem=Echem+Ener_spec_Chem(n)
     ENDDO
     WRITE(uout(3),*) 
     WRITE(uout(3),*)
     WRITE(*,'(a,E30.16)') "Ekin_v(spectral) =",Ekin   ! this might be used for testing:
     WRITE(*,'(a,E30.16)') "Etemp_v(spectral) =",Etemp ! volume avaraged energy? 
     WRITE(*,'(a,E30.16)') "Echem_v(spectral) =",Echem ! Compare with the the squared rms values 
     !                                                 ! computed in routine wrtite_diagnostics_file
     WRITE(*,'(a,E30.16)') "Emag_v(spectral) =",Emag
  ENDIF
  ! }PH
  
  ! deallocate memory
  DEALLOCATE(Ener_spec_u,Ener_spec_B,Ener_spec_Temp,Ener_spec_Chem) ! PH

END SUBROUTINE write_vertical_spectra


! Subroutine that computes the vertical power spectrum at node id zero.
FUNCTION vertical_power_spectrum(x)
  USE defprecision_module
  USE parameter_module, ONLY: Nmax
  USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  COMPLEX(kind=kr) :: x(0:,mysx_spec:,mysy_spec:)
  REAL(kind=kr)    :: vertical_power_spectrum(0:2*Nmax-1)
  REAL(kind=kr)    :: local_power_spectrum(0:2*Nmax-1)
  INTEGER(kind=ki) :: i,j,k,ierr
  REAL(kind=kr)    :: xysum,fac
  
  DO k=0,2*Nmax-1
     xysum = 0._kr
     DO j=mysy_spec,myey_spec
        DO i=mysx_spec,myex_spec
           fac = REAL((2 - DIM(1,i)),kr) 
           xysum = xysum + x(k,i,j) * CONJG(x(k,i,j)) * fac
        ENDDO
     ENDDO
     local_power_spectrum(k)=xysum
  ENDDO

  CALL MPI_REDUCE(local_power_spectrum,vertical_power_spectrum,2*Nmax,PM_MPI_FLOAT_TYPE, &
       &          MPI_SUM,0,MPI_COMM_WORLD,ierr)

END FUNCTION vertical_power_spectrum

  
