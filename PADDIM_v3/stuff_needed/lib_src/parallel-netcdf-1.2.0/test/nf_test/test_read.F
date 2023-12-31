C*********************************************************************
C   Copyright 1996, UCAR/Unidata
C   See netcdf/COPYRIGHT file for copying and redistribution conditions.
C   $Id: test_read.F 754 2009-12-30 21:19:42Z wkliao $
C*********************************************************************

C Test nfmpi_strerror.
C    Try on a bad error status.
C    Test for each defined error status.
C
        subroutine test_nfmpi_strerror()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer         number_of_messages
        parameter       (number_of_messages = 27)

        integer         i
        integer         status(number_of_messages)
        character*80    message
        character*80    msg(number_of_messages)

        data    status(1)  / NF_NOERR/
        data    status(2)  / NF_EBADID /
        data    status(3)  / NF_EEXIST /
        data    status(4)  / NF_EINVAL /
        data    status(5)  / NF_EPERM /
        data    status(6)  / NF_ENOTINDEFINE /
        data    status(7)  / NF_EINDEFINE /
        data    status(8)  / NF_EINVALCOORDS /
        data    status(9)  / NF_EMAXDIMS /
        data    status(10) / NF_ENAMEINUSE /
        data    status(11) / NF_ENOTATT /
        data    status(12) / NF_EMAXATTS /
        data    status(13) / NF_EBADTYPE /
        data    status(14) / NF_EBADDIM /
        data    status(15) / NF_EUNLIMPOS /
        data    status(16) / NF_EMAXVARS /
        data    status(17) / NF_ENOTVAR /
        data    status(18) / NF_EGLOBAL /
        data    status(19) / NF_ENOTNC /
        data    status(20) / NF_ESTS /
        data    status(21) / NF_EMAXNAME /
        data    status(22) / NF_EUNLIMIT /
        data    status(23) / NF_ENORECVARS /
        data    status(24) / NF_ECHAR /
        data    status(25) / NF_EEDGE /
        data    status(26) / NF_ESTRIDE /
        data    status(27) / NF_EBADNAME /

        data    msg(1)  / 'No error' /
        data    msg(2)  / 'Not a netCDF id' /
        data    msg(3)  / 'netCDF file exists && NC_NOCLOBBER' /
        data    msg(4)  / 'Invalid argument' /
        data    msg(5)  / 'Write to read only' /
        data    msg(6)  / 'Operation not allowed in data mode' /
        data    msg(7)  / 'Operation not allowed in define mode' /
        data    msg(8)  / 'Index exceeds dimension bound' /
        data    msg(9)  / 'NC_MAX_DIMS exceeded' /
        data    msg(10) / 'String match to name in use' /
        data    msg(11) / 'Attribute not found' /
        data    msg(12) / 'NC_MAX_ATTRS exceeded' /
        data    msg(13)
     +      / 'Not a netCDF data type or _FillValue type mismatch' /
        data    msg(14) / 'Invalid dimension id or name' /
        data    msg(15) / 'NC_UNLIMITED in the wrong index' /
        data    msg(16) / 'NC_MAX_VARS exceeded' /
        data    msg(17) / 'Variable not found' /
        data    msg(18) / 'Action prohibited on NC_GLOBAL varid' /
        data    msg(19) / 'Not a netCDF file' /
        data    msg(20) / 'In Fortran, string too short' /
        data    msg(21) / 'NC_MAX_NAME exceeded' /
        data    msg(22) / 'NC_UNLIMITED size already in use' /
        data    msg(23) / 'nc_rec op when there are no record vars' /
        data    msg(24) / 'Attempt to convert between text & numbers' /
        data    msg(25) / 'Edge+start exceeds dimension bound' /
        data    msg(26) / 'Illegal stride' /
        data    msg(27) 
     +      / 'Attribute or variable name contains illegal characters' /

C       /* Try on a bad error status */
        message = nfmpi_strerror(-666)!/* should fail */
        if (message .ne. 'Unknown Error')
     +      call errorc('nfmpi_strerror on bad error status returned: ',
     +          message)

C       /* Try on each legitimate error status */
        do 1, i=1, number_of_messages
            message = nfmpi_strerror(status(i))
            if (message .ne. msg(i))
     +          call error('nfmpi_strerror() should return "' 
     +                     // msg(i) // '"')
1       continue
        end


C Test nfmpi_open.
C If in read-only section of tests,
C    Try to open a non-existent netCDF file, check error return.
C    Open a file that is not a netCDF file, check error return.
C    Open a netCDF file with a bad mode argument, check error return.
C    Open a netCDF file with NFMPI_NOWRITE mode, try to write, check error.
C    Try to open a netcdf twice, check whether returned netcdf ids different.
C If in writable section of tests,
C    Open a netCDF file with NFMPI_WRITE mode, write something, close it.
C On exit, any open netCDF files are closed.
        subroutine test_nfmpi_open()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer err
        integer ncid
        integer ncid2
        
C       /* Try to open a nonexistent file */
        err = nfmpi_open(comm, 'tooth-fairy.nc', NF_NOWRITE,
     +                   MPI_INFO_NULL, ncid)!/* should fail */
        if (err .eq. NF_NOERR)
     +      call error
     +      ('nfmpi_open of nonexistent file should have failed')
        if (.not. NFMPI_ISSYSERR(err))
     +      call error(
     +  'nfmpi_open of nonexistent file should 
     +   have returned system error')

C       /* Open a file that is not a netCDF file. */
        err = nfmpi_open(comm, 'test_get.F', NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)!/* should fail */
        if (err .ne. NF_ENOTNC .and. err .ne. NF_EOFILE)
     +      call errore('nfmpi_open of non-netCDF file: ', err)

C       /* Open a netCDF file in read-only mode, check that write fails */
        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_redef(ncid)    !/* should fail */
        if (err .ne. NF_EPERM)
     +      call error('nfmpi_redef of read-only file should fail')
C       /* Opened OK, see if can open again and get a different netCDF ID */
        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid2)
        if (err .ne. 0) then
            call errore('nfmpi_open: ', err)
        else
            err = nfmpi_close(ncid2)
        end if
        if (ncid2 .eq. ncid)
     +      call error(
     +  'netCDF IDs for first and second
     +   nfmpi_open calls should differ')

        if (.not. readonly) then        !/* tests using netCDF scratch file */
            err = nfmpi_create(comm, scratch, NF_NOCLOBBER,
     +                         MPI_INFO_NULL, ncid2)
            if (err .ne. 0) then
                call errore('nfmpi_create: ', err)
            else 
                err = nfmpi_close(ncid2)
            end if
            err = nfmpi_open(comm, scratch, NF_WRITE, MPI_INFO_NULL,
     +                       ncid2)
            if (err .ne. 0) then
                call errore('nfmpi_open: ', err)
            else 
                err = nfmpi_close(ncid2)
            end if
            err = nfmpi_delete(scratch, MPI_INFO_NULL)
            if (err .ne. 0) 
     +          call errorc('delete of scratch file failed: ', scratch)
        end if

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


C 
C Test nfmpi_close.
C    Try to close a netCDF file twice, check whether second close fails.
C    Try on bad handle, check error return.
C    Try in define mode and data mode.
C
        subroutine test_nfmpi_close()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer err

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

C       /* Close a netCDF file twice, second time should fail */
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close failed: ', err)
        err = nfmpi_close(ncid)
        if (err .ne. NF_EBADID)
     +      call error('nfmpi_close of closed file should have failed')
        
C       /* Try with a bad netCDF ID */
        err = nfmpi_close(BAD_ID)!/* should fail */
        if (err .ne. NF_EBADID)
     +      call errore(
     +         'nfmpi_close with bad netCDF ID returned wrong error: ', 
     +         err)

C       /* Close in data mode */
        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close in data mode failed: ', err)

        if (.not. readonly) then        !/* tests using netCDF scratch file */
            err = nfmpi_create(comm, scratch, NF_NOCLOBBER,
     +                         MPI_INFO_NULL, ncid)
            if (err .ne. 0) 
     +          call errore('nfmpi_create: ', err)
            err = nfmpi_close(ncid)
            if (err .ne. 0)
     +          call errore('nfmpi_close in define mode: ', err)
            err = nfmpi_delete(scratch, MPI_INFO_NULL)
            if (err .ne. 0)
     +          call errorc('delete of scratch file failed: ', 
     +              scratch)
        end if
        end


C Test nfmpi_inq.
C    Try on bad handle, check error return.
C    Try in data mode, check returned values.
C    Try asking for subsets of info.
C If in writable section of tests,
C    Try in define mode, after adding an unlimited dimension, variable.
C On exit, any open netCDF files are closed.
        subroutine test_nfmpi_inq()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer ncid2                   !/* for scratch netCDF dataset */
        integer ndims                   !/* number of dimensions */
        integer nvars                   !/* number of variables */
        integer ngatts                  !/* number of global attributes */
        integer recdim                  !/* id of unlimited dimension */
        integer err
        integer ndims0
        integer nvars0
        integer ngatts0
        integer recdim0
        integer did
        integer vid
        integer(kind=MPI_OFFSET_KIND) length
        integer VDIMS(1)

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        
C       /* Try on bad handle */
        err = nfmpi_inq(BAD_ID, ndims, nvars, ngatts, recdim)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        
        err = nfmpi_inq(ncid, ndims, nvars, ngatts, recdim)
        if (err .ne. 0) then
            call errore('nfmpi_inq: ', err)
        else if (ndims .ne. NDIMS) then
            call errori
     +      ('nfmpi_inq: wrong number of dimensions returned: ', ndims)
        else if (nvars .ne. NVARS) then
            call errori
     +      ('nfmpi_inq: wrong number of variables returned: ', nvars)
        else if (ngatts .ne. NGATTS) then
            call errori(
     +          'nfmpi_inq: wrong number of global atts returned: ',
     +          ngatts)
        else if (recdim .ne. RECDIM) then
            call errori
     +      ('nfmpi_inq: wrong record dimension ID returned: ', recdim)
        end if

        if (.not. readonly) then        !/* tests using netCDF scratch file */
            err = nfmpi_create(comm, scratch, NF_NOCLOBBER,
     +                         MPI_INFO_NULL, ncid2)
            if (err .ne. 0) then
                call errore('nfmpi_create: ', err)
            else                !/* add dim, var, gatt, check inq */
                err = nfmpi_enddef(ncid2) !/* enter data mode */
                err = nfmpi_inq(ncid2, ndims0, nvars0, 
     +              ngatts0, recdim0)
                if (err .ne. 0)
     +              call errore('nfmpi_inq: ', err)
                err = nfmpi_redef(ncid2) !/* enter define mode */
C               /* Check that inquire still works in define mode */
                err = nfmpi_inq(ncid2, ndims, nvars, ngatts, recdim)
                if (err .ne. 0) then
                    call errore('nfmpi_inq in define mode: ', err)
                else if (ndims .ne. ndims0) then
                    call errori
     +              ('nfmpi_inq in define mode: ndims wrong, ', ndims)
                else if (nvars .ne. nvars0) then
                    call errori
     +              ('nfmpi_inq in define mode: nvars wrong, ', nvars)
                else if (ngatts .ne. ngatts0) then
                    call errori(
     +              'nfmpi_inq in define mode: ngatts wrong, ', ngatts)
                    print *, ' expected ', ngatts0
                else if (recdim .ne. recdim0) then
                    call errori
     +              ('nfmpi_inq in define mode: recdim wrong, ', recdim)
                end if

C               /* Add dim, var, global att */
                length = 1
                err = nfmpi_def_dim(ncid2, 'inqd', length, did)
                if (err .ne. 0)
     +              call errore('nfmpi_def_dim: ', err)
                VDIMS = 0
                err = nfmpi_def_var(ncid2, 'inqv', NF_FLOAT, 
     +                0, VDIMS, vid)
                if (err .ne. 0)
     +              call errore('nfmpi_def_var: ', err)

                length = len('stuff')
                err = nfmpi_put_att_text(ncid2, NF_GLOBAL, 'inqa', 
     +                                length, 'stuff')
                if (err .ne. 0)
     +              call errore('nfmpi_put_att_text: ', err)

C               /* Make sure nfmpi_inq sees the additions while in define mode */
                err = nfmpi_inq(ncid2, ndims, nvars, ngatts, recdim)
                if (err .ne. 0) then
                    call errore('nfmpi_inq in define mode: ', err)
                else if (ndims .ne. ndims0 + 1) then
                    call errori
     +              ('nfmpi_inq in define mode: ndims wrong, ', ndims)
                else if (nvars .ne. nvars0 + 1) then
                    call errori
     +              ('nfmpi_inq in define mode: nvars wrong, ', nvars)
                else if (ngatts .ne. ngatts0 + 1) then
                    call errori
     +              ('nfmpi_inq in define mode: ngatts wrong, ', ngatts)
                    print *, ' expected (added attr)', ngatts0 + 1
                end if
                err = nfmpi_enddef(ncid2)
                if (err .ne. 0)
     +              call errore('nfmpi_enddef: ', err)

C               /* Make sure nfmpi_inq stills sees additions in data mode */
                err = nfmpi_inq(ncid2, ndims, nvars, ngatts, recdim)
                if (err .ne. 0) then
                    call errore('nfmpi_inq failed in data mode: ',err)
                else if (ndims .ne. ndims0 + 1) then
                    call errori
     +              ('nfmpi_inq in define mode: ndims wrong, ', ndims)
                else if (nvars .ne. nvars0 + 1) then
                    call errori
     +              ('nfmpi_inq in define mode: nvars wrong, ', nvars)
                else if (ngatts .ne. ngatts0 + 1) then
                    call errori
     +              ('nfmpi_inq in define mode: ngatts wrong, ', ngatts)
                end if
                err = nfmpi_close(ncid2)
                err = nfmpi_delete(scratch, MPI_INFO_NULL)
                if (err .ne. 0)
     +              call errorc('delete of scratch file failed: ', 
     +                  scratch)
            end if
        end if

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_natts()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer ngatts                  !/* number of global attributes */
        integer err

        err = nfmpi_inq_natts(BAD_ID, ngatts)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_inq_natts(ncid, ngatts)
        if (err .ne. 0) then
            call errore('nfmpi_inq_natts: ', err)
        else if (ngatts .ne. NGATTS) then
            call errori
     +      ('nfmpi_inq_natts: wrong number of global atts returned, ',
     +        ngatts)
        end if
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_ndims()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer ndims
        integer err

        err = nfmpi_inq_ndims(BAD_ID, ndims)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_inq_ndims(ncid, ndims)
        if (err .ne. 0) then
            call errore('nfmpi_inq_ndims: ', err)
        else if (ndims .ne. NDIMS) then
            call errori
     +      ('nfmpi_inq_ndims: wrong number returned, ', ndims)
        end if
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_nvars()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer nvars
        integer err

        err = nfmpi_inq_nvars(BAD_ID, nvars)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_inq_nvars(ncid, nvars)
        if (err .ne. 0) then
            call errore('nfmpi_inq_nvars: ', err)
        else if (nvars .ne. NVARS) then
            call errori
     +      ('nfmpi_inq_nvars: wrong number returned, ', nvars)
        end if
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_unlimdim()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer unlimdim
        integer err

        err = nfmpi_inq_unlimdim(BAD_ID, unlimdim)
        if (err .ne. NF_EBADID)
     +      call errore('bad ncid: ', err)
        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_inq_unlimdim(ncid, unlimdim)
        if (err .ne. 0) then
            call errore('nfmpi_inq_unlimdim: ', err)
        else if (unlimdim .ne. RECDIM) then
            call errori
     +      ('nfmpi_inq_unlimdim: wrong number returned, ', unlimdim)
            print *, 'expected ', RECDIM
        end if
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_dimid()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer dimid
        integer i
        integer err

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        err = nfmpi_inq_dimid(ncid, 'noSuch', dimid)
        if (err .ne. NF_EBADDIM)
     +      call errore('bad dim name: ', err)
        do 1, i = 1, NDIMS
            err = nfmpi_inq_dimid(BAD_ID, dim_name(i), dimid)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            err = nfmpi_inq_dimid(ncid, dim_name(i), dimid)
            if (err .ne. 0) then
                call errore('nfmpi_inq_dimid: ', err)
            else if (dimid .ne. i) then
                call errori('expected ', i)
                call errori('got ', dimid)
            end if
1       continue
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_dim()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        character*(NF_MAX_NAME) name
        integer(kind=MPI_OFFSET_KIND) length
        integer intlen

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NDIMS
            err = nfmpi_inq_dim(BAD_ID, i, name, length)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            err = nfmpi_inq_dim(ncid, BAD_DIMID, name, length)
            if (err .ne. NF_EBADDIM)
     +          call errore('bad dimid: ', err)
            err = nfmpi_inq_dim(ncid, i, name, length)
            if (err .ne. 0) then
                call errore('nfmpi_inq_dim: ', err)
            else if (dim_name(i) .ne. name)  then
                call errorc('name unexpected: ', name)
                print *, ' expected ', dim_name(i),' for the ',i,'entry'
            else if (dim_len(i) .ne. length) then
                intlen = length
                call errori('size unexpected: ', intlen)
            end if
1       continue
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_dimlen()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        integer(kind=MPI_OFFSET_KIND) length
        integer intlen

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL,
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NDIMS
            err = nfmpi_inq_dimlen(BAD_ID, i, length)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            err = nfmpi_inq_dimlen(ncid, BAD_DIMID, length)
            if (err .ne. NF_EBADDIM)
     +          call errore('bad dimid: ', err)
            err = nfmpi_inq_dimlen(ncid, i, length)
            if (err .ne. 0) then
                call errore('nfmpi_inq_dimlen: ', err)
            else if (dim_len(i) .ne. length) then
                intlen = length
                call errori('size unexpected: ', intlen)
                print *, 'expected ', dim_len(i),' for the ',i,'entry'
            end if
1       continue
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_dimname()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        character*(NF_MAX_NAME)  name

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NDIMS
            err = nfmpi_inq_dimname(BAD_ID, i, name)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            err = nfmpi_inq_dimname(ncid, BAD_DIMID, name)
            if (err .ne. NF_EBADDIM)
     +          call errore('bad dimid: ', err)
            err = nfmpi_inq_dimname(ncid, i, name)
            if (err .ne. 0) then
                call errore('nfmpi_inq_dimname: ', err)
            else if (dim_name(i) .ne. name)  then
                call errorc('name unexpected: ', name)
            end if
1       continue
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_varid()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer vid
        integer i
        integer err

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

        err = nfmpi_inq_varid(ncid, 'noSuch', vid)
        if (err .ne. NF_ENOTVAR)
     +      call errore('bad ncid: ', err)

        do 1, i = 1, NVARS
            err = nfmpi_inq_varid(BAD_ID, var_name(i), vid)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            err = nfmpi_inq_varid(ncid, var_name(i), vid)
            if (err .ne. 0) then
                call errore('nfmpi_inq_varid: ', err)
            else if (vid .ne. i) then
                call errori('varid unexpected: ', vid)
            endif
1       continue

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_var()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        character*(NF_MAX_NAME) name
        integer datatype
        integer ndims
        integer dimids(MAX_RANK)
        integer na

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NVARS
            err = nfmpi_inq_var(BAD_ID, i, name, datatype, ndims, 
     +            dimids, na)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            err = nfmpi_inq_var(ncid,BAD_VARID,name,datatype,ndims,
     +                       dimids,na)
            if (err .ne. NF_ENOTVAR)
     +          call errore('bad var id: ', err)
            err = nfmpi_inq_var(ncid, i, name, datatype, ndims, dimids, 
     +                       na)
            if (err .ne. 0) then
                call errore('nfmpi_inq_var: ', err)
            else if (var_name(i) .ne. name)  then
                call errorc('name unexpected: ', name)
            else if (var_type(i) .ne. datatype) then
                call errori('type unexpected: ', datatype)
            else if (var_rank(i) .ne. ndims) then
                call errori('ndims expected: ', ndims)
            else if (.not.int_vec_eq(var_dimid(1,i),dimids,ndims)) then
                call error('unexpected dimid')
            else if (var_natts(i) .ne. na) then
                call errori('natts unexpected: ', na)
            end if
1       continue
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_vardimid()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        integer dimids(MAX_RANK)

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NVARS
            err = nfmpi_inq_vardimid(BAD_ID, i, dimids)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            err = nfmpi_inq_vardimid(ncid, BAD_VARID, dimids)
            if (err .ne. NF_ENOTVAR)
     +          call errore('bad var id: ', err)
            err = nfmpi_inq_vardimid(ncid, i, dimids)
            if (err .ne. 0) then
                call errore('nfmpi_inq_vardimid: ', err)
            else if (.not.int_vec_eq(var_dimid(1,i), dimids, 
     +               var_rank(i))) then
                call error('unexpected dimid')
                print *, ' for variable ', i
            end if
1       continue
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_varname()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        character*(NF_MAX_NAME) name

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NVARS
            err = nfmpi_inq_varname(BAD_ID, i, name)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            err = nfmpi_inq_varname(ncid, BAD_VARID, name)
            if (err .ne. NF_ENOTVAR)
     +          call errore('bad var id: ', err)
            err = nfmpi_inq_varname(ncid, i, name)
            if (err .ne. 0) then
                call errore('nfmpi_inq_varname: ', err)
            else if (var_name(i) .ne. name)  then
                call errorc('name unexpected: ', name)
            end if
1       continue
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_varnatts()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        integer na

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 0, NVARS    ! start with global attributes
            err = nfmpi_inq_varnatts(BAD_ID, i, na)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            err = nfmpi_inq_varnatts(ncid, BAD_VARID, na)
            if (err .ne. NF_ENOTVAR)
     +          call errore('bad var id: ', err)
            err = nfmpi_inq_varnatts(ncid, VARID(i), na)
            if (err .ne. 0) then
                call errore('nfmpi_inq_varnatts: ', err)
            else if (NATTS(i) .ne. na) then	! works for global attributes
                call errori('natts unexpected: ', na)
            end if
1       continue
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_varndims()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        integer ndims

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NVARS
            err = nfmpi_inq_varndims(BAD_ID, i, ndims)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            err = nfmpi_inq_varndims(ncid, BAD_VARID, ndims)
            if (err .ne. NF_ENOTVAR)
     +          call errore('bad var id: ', err)
            err = nfmpi_inq_varndims(ncid, i, ndims)
            if (err .ne. 0) then
                call errore('nfmpi_inq_varndims: ', err)
            else if (var_rank(i) .ne. ndims) then
                call errori('ndims unexpected: ', ndims)
            end if
1       continue
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_vartype()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer err
        integer datatype

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)
        do 1, i = 1, NVARS
            err = nfmpi_inq_vartype(BAD_ID, i, datatype)
            if (err .ne. NF_EBADID)
     +          call errore('bad ncid: ', err)
            err = nfmpi_inq_vartype(ncid, BAD_VARID, datatype)
            if (err .ne. NF_ENOTVAR)
     +          call errore('bad var id: ', err)
            err = nfmpi_inq_vartype(ncid, i, datatype)
            if (err .ne. 0) then
                call errore('nfmpi_inq_vartype: ', err)
            else if (var_type(i) .ne. datatype) then
                call errori('type unexpected: ', datatype)
            end if
1       continue
        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_att()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer j
        integer err
        integer t
        integer(kind=MPI_OFFSET_KIND) n

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0) 
     +      call errore('nfmpi_open: ', err)

        do 1, i = 0, NVARS
            do 2, j = 1, NATTS(i)
                err = nfmpi_inq_att(BAD_ID, i, ATT_NAME(j,i), t, n)
                if (err .ne. NF_EBADID) 
     +              call errore('bad ncid: ', err)
                err = nfmpi_inq_att
     +                     (ncid, BAD_VARID, ATT_NAME(j,i), t, n)
                if (err .ne. NF_ENOTVAR) 
     +              call errore('bad var id: ', err)
                err = nfmpi_inq_att(ncid, i, 'noSuch', t, n)
                if (err .ne. NF_ENOTATT) 
     +              call errore('Bad attribute name: ', err)
                err = nfmpi_inq_att(ncid, i, ATT_NAME(j,i), t, n)
                if (err .ne. 0) then
                    call error(nfmpi_strerror(err))
                else
                    if (t .ne. ATT_TYPE(j,i))
     +                  call error('type not that expected')
                    if (n .ne. ATT_LEN(j,i)) 
     +                  call error('length not that expected')
                end if
2           continue
1       continue

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_attlen()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer j
        integer err
        integer(kind=MPI_OFFSET_KIND) len

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

        do 1, i = 0, NVARS
            err = nfmpi_inq_attlen(ncid, i, 'noSuch', len)
            if (err .ne. NF_ENOTATT)
     +          call errore('Bad attribute name: ', err)
            do 2, j = 1, NATTS(i)
                err = nfmpi_inq_attlen(BAD_ID, i, 
     +             ATT_NAME(j,i), len)
                if (err .ne. NF_EBADID)
     +              call errore('bad ncid: ', err)
                err = nfmpi_inq_attlen
     +                (ncid, BAD_VARID, ATT_NAME(j,i), len)
                if (err .ne. NF_ENOTVAR)
     +              call errore('bad varid: ', err)
                err = nfmpi_inq_attlen(ncid, i, ATT_NAME(j,i), len)
                if (err .ne. 0) then
                    call error(nfmpi_strerror(err))
                else
                    if (len .ne. ATT_LEN(j,i))
     +                  call error('len not that expected')
                end if
2           continue
1       continue

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_atttype()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer j
        integer err
        integer datatype

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

        do 1, i = 0, NVARS
            err = nfmpi_inq_atttype(ncid, i, 'noSuch', datatype)
            if (err .ne. NF_ENOTATT)
     +          call errore('Bad attribute name: ', err)
            do 2, j = 1, NATTS(i)
                err = nfmpi_inq_atttype
     +                (BAD_ID, i, ATT_NAME(j,i), datatype)
                if (err .ne. NF_EBADID)
     +              call errore('bad ncid: ', err)
                err = nfmpi_inq_atttype(ncid, BAD_VARID, ATT_NAME(j,i), 
     +                               datatype)
                if (err .ne. NF_ENOTVAR)
     +              call errore('bad varid: ', err)
                err = nfmpi_inq_atttype
     +                (ncid, i, ATT_NAME(j,i), datatype)
                if (err .ne. 0) then
                    call error(nfmpi_strerror(err))
                else
                    if (datatype .ne. ATT_TYPE(j,i))
     +                  call error('type not that expected')
                end if
2           continue
1       continue

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_attname()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer j
        integer err
        character*(NF_MAX_NAME) name

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

        do 1, i = 0, NVARS
            err = nfmpi_inq_attname(ncid, i, BAD_ATTNUM, name)
            if (err .ne. NF_ENOTATT)
     +          call errore('Bad attribute number: ', err)
            err = nfmpi_inq_attname(ncid, i, NATTS(i)+1, name)
            if (err .ne. NF_ENOTATT)
     +          call errore('Bad attribute number: ', err)
            do 2, j = 1, NATTS(i)
                err = nfmpi_inq_attname(BAD_ID, i, j, name)
                if (err .ne. NF_EBADID)
     +              call errore('bad ncid: ', err)
                err = nfmpi_inq_attname(ncid, BAD_VARID, j, name)
                if (err .ne. NF_ENOTVAR)
     +              call errore('bad var id: ', err)
                err = nfmpi_inq_attname(ncid, i, j, name)
                if (err .ne. 0) then
                    call error(nfmpi_strerror(err))
                else
                    if (ATT_NAME(j,i) .ne. name)
     +                  call error('name not that expected')
                end if
2           continue
1       continue

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end


        subroutine test_nfmpi_inq_attid()
        use pnetcdf
        implicit        none
#include "tests.inc"
        integer ncid
        integer i
        integer j
        integer err
        integer attnum

        err = nfmpi_open(comm, testfile, NF_NOWRITE, MPI_INFO_NULL, 
     +                   ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_open: ', err)

        do 1, i = 0, NVARS
            err = nfmpi_inq_attid(ncid, i, 'noSuch', attnum)
            if (err .ne. NF_ENOTATT)
     +          call errore('Bad attribute name: ', err)
            do 2, j = 1, NATTS(i)
                err = nfmpi_inq_attid(BAD_ID, i, 
     +                ATT_NAME(j,i), attnum)
                if (err .ne. NF_EBADID)
     +              call errore('bad ncid: ', err)
                err = nfmpi_inq_attid(ncid, BAD_VARID, ATT_NAME(j,i), 
     +                             attnum)
                if (err .ne. NF_ENOTVAR)
     +              call errore('bad varid: ', err)
                err = nfmpi_inq_attid(ncid, i, 
     +               ATT_NAME(j,i), attnum)
                if (err .ne. 0) then
                    call error(nfmpi_strerror(err))
                else
                    if (attnum .ne. j)
     +                  call error('attnum not that expected')
                end if
2           continue
1       continue

        err = nfmpi_close(ncid)
        if (err .ne. 0)
     +      call errore('nfmpi_close: ', err)
        end
