! +-----------------------------------------------------------------------------+
! | The folowing subroutine writes the fields (Temp,Chem,u,curl u, B)           |
! | along the line(x=0,y=0,z) to a disk file.                                   |
! +-----------------------------------------------------------------------------+
SUBROUTINE write_z_profile(u,Temp,Chem,B,istep,t)
  USE defprecision_module
  USE defs_2D_3D_module
  USE state_module
  USE message_passing_module, ONLY : myid
  USE parameter_module, ONLY :Nmax,Gammaz,Nx,Ny,Nz, D_mag,pi,C_Lorentz
  USE mpi_transf_module, ONLY : mysx_phys,mysy_phys,mysz_phys,myez_phys,         &
  &                             mysx_spec,myex_spec,mysy_spec,myey_spec,         &
  &                             collect_xy_sum_phys
  USE diagnostics_module, ONLY : compute_z_autocorel
  IMPLICIT NONE
  TYPE(velocity)   :: u
  TYPE(field)      :: B ! PH
  TYPE(buoyancy)   :: Temp,Chem
  INTEGER(kind=ki) :: istep
  REAL(kind=kr)    :: t
  INTEGER(kind=ki) :: nprofile,nprofile_mean,ncount
  REAL(kind=kr),ALLOCATABLE :: profile(:,:) 
  REAL(kind=kr),ALLOCATABLE :: profile_mean(:,:)
  REAL(kind=kr),ALLOCATABLE :: local_z(:,:), analytical_b(:),analytical_u(:)
  REAL(kind=kr),ALLOCATABLE :: autocorel(:,:)
  COMPLEX(kind=kr),ALLOCATABLE :: temp_spec(:,:,:,:)

  REAL(kind=kr)    :: dz, zc, B_x0
  INTEGER(kind=ki) :: k

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
  
  ALLOCATE(profile(nprofile,0:Nz-1))
  ALLOCATE(profile_mean(nprofile_mean,0:Nz-1))
  ALLOCATE(local_z(nprofile,mysz_phys:myez_phys))
  ALLOCATE(autocorel(0:Nz-1,vec_x:vec_z))
  ALLOCATE(temp_spec(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,vec_x:vec_z))

  ! copy the fields(x=0,y=0,z) to local_z if you have those fields in your
  ! local memory
  IF (mysx_phys.EQ.0 .AND. mysy_phys.EQ.0) THEN
     
     local_z(1,:)=u%phys(0,0,:,vec_x)
#ifdef TWO_DIMENSIONAL
     local_z(2,:)=0._kr
#else
     local_z(2,:)=u%phys(0,0,:,vec_y)
#endif
     local_z(3,:)=u%phys(0,0,:,vec_z)
#ifdef TWO_DIMENSIONAL
     local_z(4,:)=0
     local_z(6,:)=0
#else 
     local_z(4,:)=u%curl(0,0,:,curl_x)
     local_z(6,:)=u%curl(0,0,:,curl_z)
#endif
     local_z(5,:)=u%curl(0,0,:,curl_y)
     
     ncount = 6
     
#if defined(TEMPERATURE_FIELD)
     local_z(ncount+1,:)=Temp%phys(0,0,:)
     ncount=ncount+1
#endif
#if defined(CHEMICAL_FIELD)
     local_z(ncount+1,:)=Chem%phys(0,0,:)
     ncount=ncount+1
#endif
#ifdef MAGNETIC
     local_z(ncount+1,:)=B%phys(0,0,:,vec_x)
     ncount = ncount + 1
#ifdef TWO_DIMENSIONAL
     local_z(ncount+1,:)=0._kr
     ncount = ncount+1
#else
     local_z(ncount+1,:)=B%phys(0,0,:,vec_y)
     ncount = ncount + 1
#endif
     local_z(ncount+1,:)=B%phys(0,0,:,vec_z)
     ncount = ncount + 1
#endif

  if(ncount.ne.nprofile) stop 'ncount not nprofile in write_z_prof'

  ENDIF
  
  ! collect vertical profile at x=y=0 at node 0
  CALL collect_z_profile(local_z,profile,nprofile)
  ! collect horizontal means

  CALL collect_xy_sum_phys(u%phys(:,:,:,vec_x),profile_mean(1,:),0)  
 
  CALL collect_xy_sum_phys(u%phys(:,:,:,vec_z),profile_mean(3,:),0)
  CALL collect_xy_sum_phys(u%curl(:,:,:,curl_y),profile_mean(5,:),0)
#ifdef TWO_DIMENSIONAL
  profile_mean(2,:) = 0._kr
  profile_mean(4,:) = 0._kr
  profile_mean(6,:) = 0._kr
#else 
  CALL collect_xy_sum_phys(u%phys(:,:,:,vec_y),profile_mean(2,:),0)
  CALL collect_xy_sum_phys(u%curl(:,:,:,curl_x),profile_mean(4,:),0)
  CALL collect_xy_sum_phys(u%curl(:,:,:,curl_z),profile_mean(6,:),0)
#endif

  ncount = 6 ! PH
#if defined(TEMPERATURE_FIELD)
  CALL collect_xy_sum_phys(Temp%phys(:,:,:),profile_mean(ncount+1,:),0)
  ncount = ncount + 1 
#endif

#if defined(CHEMICAL_FIELD)
  CALL collect_xy_sum_phys(Chem%phys(:,:,:),profile_mean(ncount+1,:),0)
  ncount = ncount + 1 
#endif

#ifdef MAGNETIC
  CALL collect_xy_sum_phys(B%phys(:,:,:,vec_x),profile_mean(ncount+1,:),0) ! PH
  ncount = ncount + 1
#ifdef TWO_DIMENSIONAL
  profile_mean(ncount+1,:) = 0._kr ! PH
  ncount = ncount + 1
#else
  CALL collect_xy_sum_phys(B%phys(:,:,:,vec_y),profile_mean(ncount+1,:),0) ! PH
  ncount = ncount + 1
#endif
  CALL collect_xy_sum_phys(B%phys(:,:,:,vec_z),profile_mean(ncount+1,:),0) ! PH
  ncount = ncount + 1
#endif

  if(ncount.ne.nprofile_mean) stop 'ncount not nprofile_mean in write_z_prof'


  profile_mean = profile_mean / (Nx * Ny)

  ! Compute autocorrelation functions of the velocity field

  temp_spec(:,:,:,:)=u%spec(:,:,:,:,ltime0)
  call compute_z_autocorel(temp_spec,autocorel)

  ! and write it to disk file

  dz = Gammaz / Nz
  IF (myid.EQ.0) THEN
     WRITE(uout(4),'(a,I12,a,E15.7)') "# Step=",istep," ,Time=",t 
     DO k=0,Nz-1
        WRITE(uout(4),'(E12.4,$)')  k*dz
        WRITE(uout(4),'(E12.4,$)') (profile(ncount,k),ncount=1,nprofile)
        WRITE(uout(4),'(E12.4,$)') (profile_mean(ncount,k),ncount=1,nprofile_mean)
        WRITE(uout(4),'(E12.4,$)') (autocorel(k,ncount),ncount=vec_x,vec_z)
        WRITE(uout(4),*)
     ENDDO
     WRITE(uout(4),*)
     WRITE(uout(4),*)
  ENDIF

  DEALLOCATE(profile,local_z,profile_mean,autocorel,temp_spec)

END SUBROUTINE write_z_profile


SUBROUTINE collect_z_profile(local_z,profile,nfields)
  USE defprecision_module
  USE message_passing_module, ONLY : myid,numtasks,nprocs2
  USE mpi_transf_module, ONLY : mysx_phys,myex_phys,mysy_phys,myey_phys,         &
  &                             mysz_phys,myez_phys,pencil_transpose_info_xy_zx, &
  &                             two_dim_transp_info
  USE parameter_module, ONLY : Nz
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  INTEGER(kind=ki) :: nfields
  REAL(kind=kr)    :: local_z(nfields,mysz_phys:myez_phys),profile(nfields,0:Nz-1)
  INTEGER(kind=ki) :: sz(0:nprocs2-1),ez(0:nprocs2-1) !sz(0:numtasks-1),ez(0:numtasks-1)
!  INTEGER(kind=kiMPI) :: recvcounts(0:numtasks-1),rdispls(0:numtasks-1)
  INTEGER(kind=kiMPI) :: recvcounts(0:nprocs2-1),rdispls(0:nprocs2-1)
  INTEGER(kind=kiMPI) :: ierr,comm

!  PRINT*,nprocs2
#ifdef TWO_DIMENSIONAL
  comm = two_dim_transp_info%communicator
#else
  comm = pencil_transpose_info_xy_zx%info_2nd_transpose%communicator
#endif
  IF (mysx_phys.EQ.0 .AND. mysy_phys.EQ.0) THEN
     CALL MPI_GATHER(mysz_phys,1,MPI_INTEGER,sz,1,MPI_INTEGER,0,           &
          &   comm,ierr)
     CALL MPI_GATHER(myez_phys,1,MPI_INTEGER,ez,1,MPI_INTEGER,0,           &
          &   comm,ierr)
     recvcounts(:) = (ez(:) - sz(:) + 1) * nfields
     rdispls(:) = sz(:) * nfields
     CALL MPI_GATHERV(local_z,(myez_phys-mysz_phys+1)*nfields,                &
          &    PM_MPI_FLOAT_TYPE,                                          &
          &    profile,recvcounts,rdispls,PM_MPI_FLOAT_TYPE,0,             &
          &    comm,ierr)
  ENDIF

END SUBROUTINE collect_z_profile
