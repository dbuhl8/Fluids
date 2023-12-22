! +-----------------------------------------------------------------------------+
! | The folowing subroutine writes the fields (Temp,Chem,u,curl u, B)           |
! | along the line(x=0,y,z=0) to a disk file.                                   |
! +-----------------------------------------------------------------------------+
SUBROUTINE write_y_profile(u,Temp,Chem,B,istep,t)
  USE defprecision_module
  USE state_module
  USE message_passing_module, ONLY : myid
  USE parameter_module, ONLY :Gammay,Nx,Ny,Nz, D_mag,pi,C_Lorentz
  USE mpi_transf_module, ONLY : mysx_phys,myex_phys, &
  &                             mysy_phys,myey_phys, & 
  &                             mysz_phys,myez_phys, &
  &                             mysx_spec,mysy_spec,collect_xz_sum_phys
  IMPLICIT NONE
  TYPE(velocity)   :: u
  TYPE(field)      :: B ! PH
  TYPE(buoyancy)   :: Temp,Chem
  INTEGER(kind=ki) :: istep
  REAL(kind=kr)    :: t
  INTEGER(kind=ki) :: nprofile,nprofile_mean,ncount
  REAL(kind=kr),ALLOCATABLE :: profile(:,:) 
  REAL(kind=kr),ALLOCATABLE :: profile_mean(:,:)
  REAL(kind=kr),ALLOCATABLE :: local_y(:,:)
  REAL(kind=kr)    :: dy, yc, B_x0
  INTEGER(kind=ki) :: j

  nprofile= 6
  nprofile_mean= 6
#if defined(TEMPERATURE_FIELD)
  nprofile = nprofile + 1
  nprofile_mean = nprofile_mean + 1 
#endif
#if defined(CHEMICAL_FIELD)
  nprofile = nprofile + 1
  nprofile_mean = nprofile_mean + 1 
#endif
#ifdef MAGNETIC
  nprofile = nprofile + 3
  nprofile_mean = nprofile_mean + 3
#endif


  ALLOCATE(profile(nprofile,0:Ny-1))
  ALLOCATE(profile_mean(nprofile_mean,0:Ny-1))
  ALLOCATE(local_y(nprofile,mysy_phys:myey_phys))

  ! copy the fields(x=0,y,z=0) to local_y if you have those fields in your
  ! local memory
  IF (mysx_phys.EQ.0 .AND. mysz_phys.EQ.0) THEN
     do j=mysy_phys,myey_phys     
     local_y(1,j)=u%phys(0,j,0,vec_x)
#ifdef TWO_DIMENSIONAL
     local_y(2,j)=0._kr
#else
     local_y(2,j)=u%phys(0,j,0,vec_y)
#endif
     local_y(3,j)=u%phys(0,j,0,vec_z)
#ifdef TWO_DIMENSIONAL
     local_y(4,j)=0
     local_y(6,j)=0
#else 
     local_y(4,j)=u%curl(0,j,0,curl_x)
     local_y(6,j)=u%curl(0,j,0,curl_z)
#endif
     local_y(5,j)=u%curl(0,j,0,curl_y)

     ncount = 6
     
#if defined(TEMPERATURE_FIELD)
     local_y(ncount+1,j)=Temp%phys(0,j,0)
     ncount=ncount+1
#endif
#if defined(CHEMICAL_FIELD)
     local_y(ncount+1,j)=Chem%phys(0,j,0)
     ncount = ncount + 1
#endif

#ifdef MAGNETIC
     local_y(ncount+1,j)=B%phys(0,j,0,vec_x)
     ncount = ncount + 1
#ifdef TWO_DIMENSIONAL
     local_y(ncount+1,j)=0._kr
     ncount = ncount+1
#else
     local_y(ncount+1,j)=B%phys(0,j,0,vec_y)
     ncount = ncount + 1
#endif
     local_y(ncount+1,j)=B%phys(0,j,0,vec_z)
     ncount = ncount + 1
#endif
     if(ncount.ne.nprofile) stop 'ncount not nprofile in write_y_prof'
     enddo



  ENDIF

 ! collect y-profile at x=z=0 at node 0
  CALL collect_y_profile(local_y,profile,nprofile)
  ! collect horizontal means

  CALL collect_xz_sum_phys(u%phys(:,:,:,vec_x),profile_mean(1,:),0)  
  CALL collect_xz_sum_phys(u%phys(:,:,:,vec_z),profile_mean(3,:),0)
  CALL collect_xz_sum_phys(u%curl(:,:,:,curl_y),profile_mean(5,:),0)
#ifdef TWO_DIMENSIONAL
  profile_mean(2,:) = 0._kr
  profile_mean(4,:) = 0._kr
  profile_mean(6,:) = 0._kr
#else 
  CALL collect_xz_sum_phys(u%phys(:,:,:,vec_y),profile_mean(2,:),0)
  CALL collect_xz_sum_phys(u%curl(:,:,:,curl_x),profile_mean(4,:),0)
  CALL collect_xz_sum_phys(u%curl(:,:,:,curl_z),profile_mean(6,:),0)
#endif

  ncount = 6

#if defined(TEMPERATURE_FIELD)
  CALL collect_xz_sum_phys(Temp%phys(:,:,:),profile_mean(ncount+1,:),0)
  ncount = ncount + 1 
#endif
#if defined(CHEMICAL_FIELD)
  CALL collect_xz_sum_phys(Chem%phys(:,:,:),profile_mean(ncount+1,:),0)
  ncount = ncount + 1 
#endif
#ifdef MAGNETIC
  CALL collect_xz_sum_phys(B%phys(:,:,:,vec_x),profile_mean(ncount+1,:),0) ! PH
  ncount = ncount + 1
#ifdef TWO_DIMENSIONAL
  profile_mean(ncount+1,:) = 0._kr ! PH
  ncount = ncount + 1
#else
  CALL collect_xz_sum_phys(B%phys(:,:,:,vec_y),profile_mean(ncount+1,:),0) ! PH
  ncount = ncount + 1
#endif
  CALL collect_xz_sum_phys(B%phys(:,:,:,vec_z),profile_mean(ncount+1,:),0) ! PH
  ncount = ncount + 1
#endif

  if(ncount.ne.nprofile_mean) stop 'ncount not nprofile_mean in write_y_prof'



  profile_mean = profile_mean / (Nx * Nz)

  ! and write it to disk file

  IF (myid.EQ.0) THEN
     WRITE(uout(6),'(a,I12,a,E15.7)') "# Step=",istep," ,Time=",t 
     DO j=0,Ny-1
        WRITE(uout(6),'(E12.4,$)') j*dy
        WRITE(uout(6),'(E12.4,$)') (profile(ncount,j),ncount=1,nprofile)
        WRITE(uout(6),'(E12.4,$)') (profile_mean(ncount,j),ncount=1,nprofile_mean)
        WRITE(uout(6),*)
     ENDDO
     WRITE(uout(6),*)
     WRITE(uout(6),*)
  ENDIF

  DEALLOCATE(profile,local_y,profile_mean)

END SUBROUTINE write_y_profile


SUBROUTINE collect_y_profile(local_y,profile,nfields)
  USE defprecision_module
  USE message_passing_module, ONLY : myid,numtasks,nprocs1
  USE mpi_transf_module, ONLY : mysx_phys,myex_phys,mysy_phys,myey_phys,         &
  &                             mysz_phys,myez_phys,pencil_transpose_info_zx_yz, &
  &                             two_dim_transp_info
  USE parameter_module, ONLY : Ny
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  INTEGER(kind=ki) :: nfields
  REAL(kind=kr)    :: local_y(nfields,mysy_phys:myey_phys),profile(nfields,0:Ny-1)
  INTEGER(kind=ki) :: sy(0:nprocs1-1),ey(0:nprocs1-1) 
!  INTEGER(kind=ki) :: sy(0:numtasks-1),ey(0:numtasks-1)
!  INTEGER(kind=kiMPI) :: recvcounts(0:numtasks-1),rdispls(0:numtasks-1)
  INTEGER(kind=kiMPI) :: recvcounts(0:nprocs1-1),rdispls(0:nprocs1-1)
  INTEGER(kind=kiMPI) :: ierr,comm

!  PRINT*,nprocs2

#ifdef TWO_DIMENSIONAL
  comm = two_dim_transp_info%communicator
#else
  comm = pencil_transpose_info_zx_yz%info_1st_transpose%communicator
#endif
  IF (mysx_phys.EQ.0 .AND. mysz_phys.EQ.0) THEN
     CALL MPI_GATHER(mysy_phys,1,MPI_INTEGER,sy,1,MPI_INTEGER,0,           &
          &   comm,ierr)
     CALL MPI_GATHER(myey_phys,1,MPI_INTEGER,ey,1,MPI_INTEGER,0,           &
          &   comm,ierr)
     recvcounts(:) = (ey(:) - sy(:) + 1) * nfields
     rdispls(:) = sy(:) * nfields
     CALL MPI_GATHERV(local_y,(myey_phys-mysy_phys+1)*nfields,                &
          &    PM_MPI_FLOAT_TYPE,                                          &
          &    profile,recvcounts,rdispls,PM_MPI_FLOAT_TYPE,0,             &
          &    comm,ierr)
  ENDIF

END SUBROUTINE collect_y_profile
