SUBROUTINE write_compressed_file(u,Temp,Chem,istep,t,dt)
  USE defprecision_module
  USE state_module, ONLY: buoyancy,velocity
  USE mpi_transf_module, ONLY : mysy_phys,myey_phys,mysz_phys,myez_phys,  &
  &                             collect_phys
  USE message_passing_module, ONLY : myid
  USE parameter_module, ONLY: Nx,Ny,Nz
  IMPLICIT NONE
  TYPE(velocity)             :: u
  TYPE(buoyancy)             :: Temp,Chem
  REAL(kind=kr)              :: t,dt
  INTEGER(kind=ki)           :: istep
  REAL(kind=krs)             :: local_sp(0:Nx-1,mysy_phys:myey_phys, &
  &                                             mysz_phys:myez_phys)
  REAL(kind=krs),ALLOCATABLE :: work(:,:,:)
  INTEGER(kind=ki)           :: ierr,jcwrite,quality
#ifdef TWO_DIMENSIONAL
  REAL(kind=krs)             :: work2(nx*nz)
  INTEGER(kind=ki)           :: i
#endif
!
  IF (write_compressed_fields) THEN 
     ! allocate temporary work space
     IF (myid.EQ.0) THEN 
        ALLOCATE(work(0:Nx-1,0:Ny-1,0:Nz-1))
     ELSE 
        ALLOCATE(work(1,1,1)) ! some compilers complain if not allocated
     ENDIF
     ! this quality parameter determines the accuracy of the compressed data
     quality=80
     ! convert data to real precision and pass to the jc-library
 
#ifndef TWO_DIMENSIONAL 
! ----------------------------------- 3D case ---------------------------------
#ifdef TEMPERATURE_FIELD
     local_sp = REAL(Temp%phys,kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0)                                                           &                                                  
     &   ierr = jcwrite(1,JPTEMP, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
     &                  work, Nx, Ny, Nz,quality)
#endif

#ifdef CHEMICAL_FIELD
     local_sp = REAL(Chem%phys,kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0)                                                           &  
     &   ierr = jcwrite(1,JPCHEM, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
     &                  work, Nx, Ny, Nz,quality)
#endif

     local_sp = REAL(u%phys(:,:,:,vec_x),kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0)                                                           &  
     &   ierr = jcwrite(1,JPVELX, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
     &                  work, Nx, Ny, Nz,quality)

     local_sp = REAL(u%phys(:,:,:,vec_y),kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0)                                                           &  
     &   ierr = jcwrite(1,JPVELY, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
     &                  work, Nx, Ny, Nz,quality)

     local_sp = REAL(u%phys(:,:,:,vec_z),kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0)                                                           &  
     &   ierr = jcwrite(1,JPVELZ, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
     &                  work, Nx, Ny, Nz,quality)

     local_sp = REAL(u%curl(:,:,:,curl_x),kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0)                                                           &  
     &   ierr = jcwrite(1,JPBX, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
     &                  work, Nx, Ny, Nz,quality)

     local_sp = REAL(u%curl(:,:,:,curl_y),kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0)                                                           &  
     &   ierr = jcwrite(1,JPBY, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
     &                  work, Nx, Ny, Nz,quality)

     local_sp = REAL(u%curl(:,:,:,curl_z),kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0)                                                           &  
     &   ierr = jcwrite(1,JPBZ, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
     &                  work, Nx, Ny, Nz,quality)

#endif

#ifdef TWO_DIMENSIONAL 
! ----------------------------------- 2D case ---------------------------------
! Use Joerg's raterfari tool for visualization
     ! Vorticity field
     local_sp = REAL(u%curl(:,:,:,curl_y),kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0) THEN
        DO i=1,nx*nz
           work2(i) = work( (i-1)/nz,0,MOD(i-1,nz))
        ENDDO
        IF (myid.EQ.0)                                                           &  
        &     ierr = jcwrite(1,JCPSI, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
        &                    work2, Nx, Nz, 1,quality)
        IF (ierr.NE.0) WRITE(*,*) "Error jcwrite"
     ENDIF
     ! Velocity field
     local_sp = REAL(u%phys(:,:,:,vec_x),kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0) THEN
        DO i=1,nx*nz
           work2(i) = work( (i-1)/nz,0,MOD(i-1,nz))
        ENDDO
        IF (myid.EQ.0)                                                           &  
        &     ierr = jcwrite(1,JCPSIDX, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
        &                    work2, Nx, Nz, 1,quality)
        IF (ierr.NE.0) WRITE(*,*) "Error jcwrite"
     ENDIF
     local_sp = REAL(u%phys(:,:,:,vec_z),kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0) THEN
        DO i=1,nx*nz
           work2(i) = work( (i-1)/nz,0,MOD(i-1,nz))
        ENDDO
        IF (myid.EQ.0)                                                           &  
        &     ierr = jcwrite(1,JCPSIDY, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
        &                    work2, Nx, Nz, 1,quality)
        IF (ierr.NE.0) WRITE(*,*) "Error jcwrite"
     ENDIF
#ifdef TEMPERATURE_FIELD
     ! Temperature field
     local_sp = REAL(Temp%phys,kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0) THEN
        DO i=1,nx*nz
           work2(i) = work( (i-1)/nz,0,MOD(i-1,nz))
        ENDDO
        IF (myid.EQ.0)                                                             &  
        &     ierr = jcwrite(1,JCTEMP, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
        &                    work2, Nx, Nz, 1,quality)
        IF (ierr.NE.0) WRITE(*,*) "Error jcwrite"
     ENDIF
#endif
#ifdef CHEMICAL_FIELD
     local_sp = REAL(Chem%phys,kind=krs)
     CALL collect_phys(local_sp,work,0)
     IF (myid.EQ.0) THEN
        DO i=1,nx*nz
           work2(i) = work( (i-1)/nz,0,MOD(i-1,nz))
        ENDDO
        IF (myid.EQ.0)                                                             &  
        &     ierr = jcwrite(1,JCCHEM, istep, REAL(t,kind=krs), REAL(dt,kind=krs), &
        &                    work2, Nx, Nz, 1,quality)
        IF (ierr.NE.0) WRITE(*,*) "Error jcwrite"
     ENDIF
#endif

#endif

     DEALLOCATE(work)
  ENDIF

END SUBROUTINE write_compressed_file
