! File: read_stochastic_forcing.f90
! Author: Dante Buhl
! Purpose: Read stochastic forcing from text file and place in a
! fortran array

subroutine read_stochastic_forcing(restarted)

    USE defprecision_module
    USE parameter_module
    USE message_passing_module, ONLY: myid
    USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec
    
    
#ifdef MPI_MODULE
    USE MPI
#else
    INCLUDE "mpif.h"
#endif

    implicit none
    
    character :: restarted
    logical :: gprestart
    integer(kind=ki) :: i, j, readNumColumns, ierr, rs, numTrashLines
    real(kind=kr) :: hkx, hky
    real(kind=kr) :: max_forced_wavenumber, timeStop
    integer(kind=ki) :: number_of_points_in_window = 10_ki
    integer(kind=ki) :: number_of_points_between_window_update = 250_ki
    real(kind=kr) :: delta_t = 0.5_kr
    real(kind=kr) :: gaussian_timescale = 1._kr
    real(kind=kr) :: relative_eigenvalue_tolerance = 10_kr**-4
    real(kind=kr) :: last_timestep_of_forcing
    character*50 :: string
    
    allocate(waveNumMap(mysx_spec:myex_spec, mysy_spec:myey_spec))

    headerFormatReal = "('# ', A, F15.8)"
    headerFormatInt = "('# ', A, I6)"
  
 
    NAMELIST /forcing_values/ max_forced_wavenumber
    NAMELIST /forcing_values/ number_of_points_in_window
    NAMELIST /forcing_values/ number_of_points_between_window_update 
    NAMELIST /forcing_values/ delta_t                                
    NAMELIST /forcing_values/ gaussian_timescale                     
    NAMELIST /forcing_values/ relative_eigenvalue_tolerance          

    IF (myid.EQ.0) THEN
    ! Read data from given namelist
    ! -----------------------------
!        OPEN(25, file='forcing_parameters')
        READ(*,NML=forcing_values,IOSTAT=ierr)
        IF (ierr .NE. 0 ) THEN
            PRINT*,myid,'READ ERROR IN INPUT NAMELIST'
            STOP
        ENDIF
!        CLOSE(25)
    
        KMAX_forcing                =  max_forced_wavenumber
        window_pts                  =  number_of_points_in_window
        window_skip                 =  number_of_points_between_window_update
        tstep                       =  delta_t
        gaussian_tmscl              =  gaussian_timescale 
        tol                         =  relative_eigenvalue_tolerance
    
    ENDIF
    ! Broadcast this to the rest of the processors. 
    CALL MPI_BCAST(KMAX_forcing          ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(window_pts            ,1,   MPI_INTEGER   ,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(window_skip           ,1,   MPI_INTEGER   ,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(gaussian_tmscl        ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(tol                   ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(tstep                 ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)

    numGPcolumns = 0_ki
    waveNumMap = 1_ki
    numGProws = NINT((Nsteps*dt_max)/tstep, ki) + 1_ki
   
    ! Storing values in waveNumMap
    ! ----------------------------------------------------------------------------------------- 
    if(Mmax .ne. 0) then
        do j=mysy_spec,myey_spec
            hky = ky(j)
            do i=mysx_spec,myex_spec
                hkx = kx(i)
                if((sqrt(hkx**2 + hky**2) .le. KMAX_forcing) .and. (hkx + hky .gt. 0_ki)) then 
                    numGPcolumns = numGPcolumns + 1_ki
                    waveNumMap(i, j) = numGPcolumns
                end if 
            end do
        end do
    else 
        do i = mysx_spec, myex_spec
            if((kx(i) .le. KMAX_forcing) .and. (kx(i) .gt. 0_ki)) then
                numGPcolumns = numGPcolumns + 1_ki
                waveNumMap(i, :) = numGPcolumns
            end if 
        end do
    end if
    print *, "cpu "//trim(str(myid))//": allocated wavenumber map"
    ! ----------------------------------------------------------------------------------------- 

    ! ----------------------------------------------------------------------------------------- 
    allocate(gpForcingVals(numGProws, numGPcolumns))
    allocate(gpTimeVals(numGProws)) 
    allocate(interpSlopeVals(numGProws-1_ki, numGPcolumns))
    dataFormat = "("//trim(str(numGPcolumns+1))//"(F7.4, '    '))"
    ! ----------------------------------------------------------------------------------------- 


    ! If this processor needs forcing, then do the forcing 
    ! ----------------------------------------------------------------------------------------- 
    if(numGPcolumns .gt. 0) then
        print *, "cpu "//trim(str(myid))//": checking if there is already a forcing file"
    
        timeStop = Nsteps*dt_max + 10.0_kr

        ! Make sure there is a forcing file for this processor with enough columns
        ! ------------------------------------------------------------------------------------- 
        inquire(file="forcing_data/forcing"//trim(str(myid))//".dat", iostat=rs)

        if((restarted .eq. "N") .or. (rs .ne. 0)) then 

            CALL gaussian(numGPcolumns, window_pts, window_skip, timeStop, tstep, &
                          gaussian_tmscl, tol, myid, numTrashLines)
        else
            print *, "cpu "//trim(str(myid))//": reading parameters from forcing file"

            open(20, file="forcing_data/forcing"//trim(str(myid))//".dat")
                read(20, headerFormatInt)  string, readNumColumns
                read(20, headerFormatInt)  string, number_of_points_in_window
                read(20, headerFormatInt)  string, number_of_points_between_window_update
                read(20, headerFormatReal) string, gaussian_timescale
                read(20, headerFormatReal) string, delta_t
                read(20, headerFormatReal) string, relative_eigenvalue_tolerance
                read(20, headerFormatReal) string, last_timestep_of_forcing
                read(20, headerFormatInt)  string, numTrashLines
            close(20)


            print *, "cpu "//trim(str(myid))//": checking if parameters match", readNumColumns, &
                           numTrashLines, delta_t, gaussian_timescale

            gprestart = .false.
            if (readNumColumns .lt. numGPcolumns) then
                gprestart = .true.
            else if (number_of_points_in_window .ne. window_pts) then
                gprestart = .true.
            else if (number_of_points_between_window_update .ne. window_skip) then
                gprestart = .true.
            else if (gaussian_timescale .ne. gaussian_tmscl) then
                gprestart = .true.
            else if (delta_t .ne. tstep) then
                gprestart = .true.
            else if (relative_eigenvalue_tolerance .ne. tol) then
                gprestart = .true.
            else if (last_timestep_of_forcing .ne. timeStop) then
                gprestart = .true.
            end if

            if(gprestart) then
                print *, "cpu "//trim(str(myid))//": regenerating gp files"
                CALL gaussian(numGPcolumns, window_pts, window_skip, timeStop, &
                              tstep, gaussian_tmscl, tol, myid, numTrashLines)
            end if

        end if
        ! ------------------------------------------------------------------------------------- 
        
        
        ! Read Values from the File
        ! ------------------------------------------------------------------------------------- 
        open(21, file="forcing_data/forcing"//trim(str(myid))//".dat")
           
            print *, "cpu "//trim(str(myid))//": read_stoch, reading forcing data, numtrashlines:", numTrashLines
            do i=1, numTrashLines
                read(21, *)
            end do
        
            print *, "cpu "//trim(str(myid))//": read_stoch, passed header lines."

            do i=1, numGProws
                read(21,dataFormat) gpTimeVals(i), gpForcingVals(i, :)
                if(myid .eq. 0) then
                    print *, "cpu "//trim(str(myid))//": read_stoch, loop,", i
                end if
            end do
            
            print *, "cpu "//trim(str(myid))//": read_stoch, passed data lines"
        close(21)
        ! ------------------------------------------------------------------------------------- 
       
        ! Compute Slope Values for Interpolation 
        ! ------------------------------------------------------------------------------------- 
        do i=1, numGProws-1
            interpSlopeVals(i, :) = (gpForcingVals(i+1,:) - gpForcingVals(i,:))/&
                                    &(gpTimeVals(i+1)-gpTimeVals(i))
        end do
        ! ------------------------------------------------------------------------------------- 


    end if
    ! ----------------------------------------------------------------------------------------- 

    print *, "cpu "//trim(str(myid))//": completed read_stochastic_forcing.f90 call"

 
    contains
    
    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

end subroutine read_stochastic_forcing
