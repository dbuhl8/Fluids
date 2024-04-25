

program shearSolve

    use LinAl

    implicit none

    integer, parameter :: M = 50, Nmax = 1, nk = 1
    integer, parameter :: LDA = 5*(2*Nmax+1)
    real, dimension(nk, 1) :: kz, kx
    real :: ky
    real :: kxmin = 0.01, kxmax = 1.5, kzmin = 0.01, kzmax = 1.5
    real, dimension(LDA) :: ALPHAR, ALPHAI, BETA
    integer :: i, j, k, indu, indv, indw, indt, indp, indm
    real :: Pe, Re, Ri, f, dispterm
    real :: uf, vf, wf, tf, pf, us, alpha
    character :: jobvl = "N", jobvr = "V"
    real, dimension(LDA, LDA) :: A, B, D, V, VL
    real, dimension(LDA, LDA) :: VR
    real, dimension(nk, nk) :: lambda
    real, allocatable :: work(:)
    integer, dimension(2, 1) :: maxidx
    integer :: lwork, info

    Pe = 1000.0
    Re = 1000.0
    Ri = 100.0

    ky = 1.0
    f = 0.0
   
!   do i = 1, nk
!       kx(i, 1) = kxmin + (i-1)*(kxmax-kxmin)/(nk-1) 
!       do j = 1, nk
!           kz(j, 1) = kzmin + (j-1)*(kzmax - kzmin)/(nk-1)

!           A = 0.0
!           B = 0.0
!   
!           do k = -Nmax, Nmax
!               indu = k + Nmax + 1
!               indv = 2*Nmax + 1 + indu
!               indw = 2*Nmax + 1 + indv
!               indt = 2*Nmax + 1 + indw
!               indp = 2*Nmax + 1 + indt

!               dispterm = kx(i, 1)**2 + ((f+k)**2)*ky**2 + kz(j, 1)**2
!               if (k > -Nmax) then
!                   A(indu, indu-1) = -kx(i, 1)/2
!                   A(indu, indv-1) = -ky/2
!                   A(indw, indw-1) = -kx(i, 1)/2
!                   A(indt, indt-1) = -kx(i, 1)/2
!                   A(indv, indv-1) = -kx(i, 1)/2
!               end if
!               if( k < Nmax) then
!                   A(indu, indu+1) = kx(i, 1)/2
!                   A(indu, indv+1) = -ky/2
!                   A(indw, indw+1) = kx(i, 1)/2
!                   A(indt, indt+1) = kx(i, 1)/2
!                   A(indv, indv+1) = kx(i, 1)/2
!               end if 

!               ! Dissipitive Terms
!               A(indu,indu)= -dispterm/Re;
!               A(indw,indw)= -dispterm/Re;
!               A(indv,indv)= -dispterm/Re;
!               A(indt,indt)= -dispterm/Pe;

!               ! Stratification Terms here
!               A(indt,indw)=-1.0;
!               A(indw,indt)= Ri;

!               ! Pressure Terms
!               A(indu,indp)=-kx(i, 1);  ! avoid complex matrix by defining variable ip
!               A(indv,indp)=-ky*(f+k);
!               A(indw,indp)=-kz(j, 1);

!               ! Eigenvalue Matrix B, 1 on diagonal
!               B(indu,indu)=1.0;
!               B(indv,indv)=1.0;
!               B(indw,indw)=1.0;
!               B(indt,indt)=1.0;
!               B(indp,indp) = 0.0;
!               
!               ! Continuity Equation
!               A(indp,indu) = kx(i, 1);
!               A(indp,indv) = ky*(f+k);
!               A(indp,indw) = kz(j, 1);

!           end do

!           info = 0
!           VL = 0.0
!           VR = 0.0
!           ALPHAR = 0.0
!           ALPHAI = 0.0
!           BETA = 0.0
!           lwork = -1
!           allocate(work(1))
!           work = 0.0
!           call dggev(jobvl, jobvr, LDA, A, LDA, B, LDA, ALPHAR, ALPHAI, BETA, VL, LDA, VR, LDA, WORK, LWORK, info)
!           lwork = work(1)
!           deallocate(work)  
!           allocate(work(lwork))
!           call dggev(jobvl, jobvr, LDA, A, LDA, B, LDA, ALPHAR, ALPHAI, BETA, VL, LDA, VR, 2*LDA, WORK, LWORK, info)
!           deallocate(work)
!           do k = 1, LDA
!               if(BETA(k) > 10.d-6) then
!                   D(k, k) = ALPHAR(k)/BETA(k)
!               end if
!           end do
!           indm = 1 
!           do k = 1, LDA
!               if(D(k, k) .le. 0) then
!                   D(k, k) = 0.;
!               end if
!               if(D(k, k) .ge. 100.0) then
!                   D(k, k) = 0.;
!               end if
!               if(D(k, k) > D(indm, indm)) then
!                   indm = k
!               end if
!           end do
!           lambda(i, j) = D(indm, indm)
!       end do
!   end do

!   

!   call printmat(lambda, nk, nk)
!   call printmat(kx, nk, 1)
!   call printmat(kz, nk, 1)

!   open(10, file="kz.dat")
!       do i = 1, nk
!          write (10, *) kz(i, 1) 
!       end do
!   close(10)
!   open(11, file="kx.dat")
!       do i = 1, nk
!          write (11, *) kx(i, 1) 
!       end do
!   close(11)
!   open(12, file="lambda.dat")
!       do i = 1, nk
!          write (12, *) lambda(i, :)
!       end do
!   close(12)
!   open(13, file="nk.dat")
!       write (13, *) nk
!   close(13)


    kx = 0.5
    kz = 0.0

    !solve the eigenproblem 1 more time to recover the eigenvectors
            A = 0.0
            B = 0.0
            i = 1
            j = 1
    
            do k = -Nmax, Nmax
                indu = k + Nmax + 1
                indv = 2*Nmax + 1 + indu
                indw = 2*Nmax + 1 + indv
                indt = 2*Nmax + 1 + indw
                indp = 2*Nmax + 1 + indt

                dispterm = kx(i, 1)**2 + ((f+k)**2)*ky**2 + kz(j, 1)**2
                if (k > -Nmax) then
                    A(indu, indu-1) = -kx(i, 1)/2
                    A(indu, indv-1) = -ky/2
                    A(indw, indw-1) = -kx(i, 1)/2
                    A(indt, indt-1) = -kx(i, 1)/2
                    A(indv, indv-1) = -kx(i, 1)/2
                end if
                if( k < Nmax) then
                    A(indu, indu+1) = kx(i, 1)/2
                    A(indu, indv+1) = -ky/2
                    A(indw, indw+1) = kx(i, 1)/2
                    A(indt, indt+1) = kx(i, 1)/2
                    A(indv, indv+1) = kx(i, 1)/2
                end if 

                ! Dissipitive Terms
                A(indu,indu)= -dispterm/Re;
                A(indw,indw)= -dispterm/Re;
                A(indv,indv)= -dispterm/Re;
                A(indt,indt)= -dispterm/Pe;

                ! Stratification Terms here
                A(indt,indw)=-1.0;
                A(indw,indt)= Ri;

                ! Pressure Terms
                A(indu,indp)=-kx(i, 1);  ! avoid complex matrix by defining variable ip
                A(indv,indp)=-ky*(f+k);
                A(indw,indp)=-kz(j, 1);

                ! Eigenvalue Matrix B, 1 on diagonal
                B(indu,indu)=1.0;
                B(indv,indv)=1.0;
                B(indw,indw)=1.0;
                B(indt,indt)=1.0;
                B(indp,indp) = 0.0;
                
                ! Continuity Equation
                A(indp,indu) = kx(i, 1);
                A(indp,indv) = ky*(f+k);
                A(indp,indw) = kz(j, 1);

            end do

            info = 0
            VL = 0.0
            VR = 0.0
            ALPHAR = 0.0
            ALPHAI = 0.0
            BETA = 0.0
            lwork = -1
            allocate(work(1))
            work = 0.0
            call dggev(jobvl, jobvr, LDA, A, LDA, B, LDA, ALPHAR, ALPHAI, BETA, VL, LDA, VR, LDA, WORK, LWORK, info)
            lwork = work(1)
            deallocate(work)  
            allocate(work(lwork))
            call dggev(jobvl, jobvr, LDA, A, LDA, B, LDA, ALPHAR, ALPHAI, BETA, VL, LDA, VR, LDA, WORK, LWORK, info)
            deallocate(work)
            do k = 1, LDA
                if(BETA(k) > 10.d-6) then
                    D(k, k) = ALPHAR(k)/BETA(k)
                end if
            end do
            indm = 1 
            do k = 1, LDA
                if(D(k, k) .le. 0) then
                    D(k, k) = 0.;
                end if
                if(D(k, k) .ge. 100.0) then
                    D(k, k) = 0.;
                end if
                if(D(k, k) > D(indm, indm)) then
                    indm = k
                end if
            end do
            lambda(i, j) = D(indm, indm)

            call printmat(vr(:, indm:indm), LDA, 1)
            

  
end program shearSolve




