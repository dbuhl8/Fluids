SUBROUTINE write_shell_spectra(u,Temp,Chem,B,istep,t)
  USE defprecision_module
  USE state_module, ONLY: buoyancy, velocity,field,ltime0 ! PH
  USE message_passing_module, ONLY : myid
  USE parameter_module, ONLY: kx, ky, kz,Nmax, alpha, beta, gamma
  USE mpi_transf_module,  ONLY:  mysx_spec,myex_spec,mysy_spec,myey_spec
  IMPLICIT NONE
  TYPE(velocity)   :: u
  TYPE(field)      :: B ! PH
  TYPE(buoyancy)   :: Temp,Chem
  INTEGER(kind=ki) :: istep
  REAL(kind=kr)    :: t, bin_sz
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

  bin_sz = min(alpha, beta, gamma)

  Ener_spec_u(:,1) = shell_power_spectrum(u%spec(:,:,:,vec_x,ltime0))
#ifdef TWO_DIMENSIONAL
  Ener_spec_u(:,2) = 0._kr
#else
  Ener_spec_u(:,2) = shell_power_spectrum(u%spec(:,:,:,vec_y,ltime0))
#endif
  Ener_spec_u(:,3) = shell_power_spectrum(u%spec(:,:,:,vec_z,ltime0))
#ifdef TEMPERATURE_FIELD
  Ener_spec_Temp   = shell_power_spectrum(Temp%spec(:,:,:,ltime0))
#else
  Ener_spec_Temp   = 0._kr
#endif
#ifdef CHEMICAL_FIELD
  Ener_spec_Chem   = shell_power_spectrum(Chem%spec(:,:,:,ltime0))
#else
  Ener_spec_Chem   = 0._kr
#endif

#ifdef MAGNETIC
  Ener_spec_B(:,1) = shell_power_spectrum(B%spec(:,:,:,vec_x,ltime0))
#ifdef TWO_DIMENSIONAL
  Ener_spec_B(:,2) = 0._kr
#else
  Ener_spec_B(:,2) = shell_power_spectrum(B%spec(:,:,:,vec_y,ltime0))
#endif
  Ener_spec_B(:,3) = shell_power_spectrum(B%spec(:,:,:,vec_z,ltime0))
#else
  Ener_spec_B = 0._kr
#endif

  ! write spectrum
  ! PH{
  IF (myid.EQ.0) THEN 
     WRITE (uout(5),'(a,I8,a,E20.7)') '# Timstep =',istep,' time =',t
     Ekin=0._kr ! used for testing: see below
     Etemp=0._kr
     EChem=0._kr
     DO n=0,2*Nmax-1
        WRITE(uout(5),'(11E20.7)') n*bin_sz,                      &
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
     WRITE(uout(5),*) 
     WRITE(uout(5),*)
     WRITE(*,'(a,E30.16)') "Ekin_v(spectral) =",Ekin   ! this might be used for testing:
     WRITE(*,'(a,E30.16)') "Etemp_v(spectral) =",Etemp ! volume avaraged energy? 
     WRITE(*,'(a,E30.16)') "Echem_v(spectral) =",Echem ! Compare with the the squared rms values 
     !                                                 ! computed in routine wrtite_diagnostics_file
     WRITE(*,'(a,E30.16)') "Emag_v(spectral) =",Emag
  ENDIF
  ! }PH
  
  ! deallocate memory
  DEALLOCATE(Ener_spec_u,Ener_spec_B,Ener_spec_Temp,Ener_spec_Chem) ! PH

END SUBROUTINE write_shell_spectra


! Subroutine that computes the shell power spectrum at node id zero.
FUNCTION shell_power_spectrum(x)
  USE defprecision_module
  USE parameter_module, ONLY: Nmax, kx, ky, kz, alpha, beta, gamma
  USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  COMPLEX(kind=kr) :: x(0:,mysx_spec:,mysy_spec:)
  REAL(kind=kr)    :: shell_power_spectrum(0:2*Nmax-1)
  REAL(kind=kr)    :: local_shell_spectrum(0:2*Nmax-1)
  INTEGER(kind=ki) :: i,j,k,ierr,Num
  REAL(kind=kr)    :: fac,bin_sz,k_net

  local_shell_spectrum = 0._kr
  bin_sz = min(alpha, beta, gamma)

  DO j=mysy_spec,myey_spec
    DO i=mysx_spec,myex_spec
      DO k=0, 2*Nmax-1
        k_net = sqrt(kx(i)*kx(i) + ky(j)*ky(j) + kz(k)*kz(k))
        Num = floor(k_net/bin_sz-0.5)+1
        fac = REAL((2 - DIM(1,i)),kr)
        IF(Num.gt.2*Nmax-1) Num=2*Nmax-1 
        local_shell_spectrum(Num) = local_shell_spectrum(Num) + x(k,i,j) * CONJG(x(k,i,j)) * fac
        ENDDO
     ENDDO
  ENDDO

  CALL MPI_REDUCE(local_shell_spectrum,shell_power_spectrum,2*Nmax,PM_MPI_FLOAT_TYPE, &
       &          MPI_SUM,0,MPI_COMM_WORLD,ierr)

END FUNCTION shell_power_spectrum

  
