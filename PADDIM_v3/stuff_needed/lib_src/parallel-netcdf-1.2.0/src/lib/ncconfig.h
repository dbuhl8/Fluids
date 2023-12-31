/* src/lib/ncconfig.h.  Generated from ncconfig.in by configure.  */
/* src/lib/ncconfig.in.  Generated from configure.in by autoheader.  */

/* Define to one of `_getb67', `GETB67', `getb67' for Cray-2 and Cray-YMP
   systems. This function is required for `alloca.c' support on those systems.
   */
/* #undef CRAY_STACKSEG_END */

/* Define to 1 if using `alloca.c'. */
/* #undef C_ALLOCA */

/* Define if able to support CDF-5 file format */
#define ENABLE_CDF5 

/* Define if able to support nonblocking routines */
#define ENABLE_NONBLOCKING 

/* Define if Fortran names are lower case */
/* #undef F77_NAME_LOWER */

/* Define if Fortran names are lower case with two trailing underscore2 */
/* #undef F77_NAME_LOWER_2USCORE */

/* Define if Fortran names are lower case with one trailing underscore */
#define F77_NAME_LOWER_USCORE 

/* Define if Fortran names are uppercase */
/* #undef F77_NAME_UPPER */

/* Define to 1 if you have `alloca', as a function or macro. */
#define HAVE_ALLOCA 1

/* Define to 1 if you have <alloca.h> and it should be used (not on Ultrix).
   */
#define HAVE_ALLOCA_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* available */
#define HAVE_MPI_CHARACTER 1

/* available */
#define HAVE_MPI_COMBINER_DARRAY 1

/* available */
#define HAVE_MPI_COMBINER_DUP 1

/* available */
#define HAVE_MPI_COMBINER_F90_COMPLEX 1

/* available */
#define HAVE_MPI_COMBINER_F90_INTEGER 1

/* available */
#define HAVE_MPI_COMBINER_F90_REAL 1

/* available */
#define HAVE_MPI_COMBINER_HINDEXED_INTEGER 1

/* available */
#define HAVE_MPI_COMBINER_HVECTOR_INTEGER 1

/* available */
#define HAVE_MPI_COMBINER_INDEXED_BLOCK 1

/* available */
#define HAVE_MPI_COMBINER_RESIZED 1

/* available */
#define HAVE_MPI_COMBINER_STRUCT_INTEGER 1

/* available */
#define HAVE_MPI_COMBINER_SUBARRAY 1

/* available */
#define HAVE_MPI_COMPLEX16 1

/* available */
#define HAVE_MPI_COMPLEX32 1

/* available */
#define HAVE_MPI_COMPLEX8 1

/* available */
#define HAVE_MPI_DOUBLE_PRECISION 1

/* Define to 1 if you have the `MPI_Info_dup' function. */
#define HAVE_MPI_INFO_DUP 1

/* Define to 1 if you have the `MPI_Info_free' function. */
#define HAVE_MPI_INFO_FREE 1

/* available */
#define HAVE_MPI_INTEGER 1

/* available */
#define HAVE_MPI_INTEGER1 1

/* available */
#define HAVE_MPI_INTEGER16 1

/* available */
#define HAVE_MPI_INTEGER2 1

/* available */
#define HAVE_MPI_INTEGER4 1

/* available */
#define HAVE_MPI_INTEGER8 1

/* available */
#define HAVE_MPI_LB 1

/* available */
#define HAVE_MPI_REAL 1

/* available */
#define HAVE_MPI_REAL16 1

/* available */
#define HAVE_MPI_REAL4 1

/* available */
#define HAVE_MPI_REAL8 1

/* Define to 1 if you have the `MPI_Request_get_status' function. */
#define HAVE_MPI_REQUEST_GET_STATUS 1

/* RESIZED combiner works */
#define HAVE_MPI_RESIZED_SUPPORT 1

/* Define to 1 if you have the `MPI_Type_dup' function. */
#define HAVE_MPI_TYPE_DUP 1

/* available */
#define HAVE_MPI_UB 1

/* Define to 1 if the system has the type `ptrdiff_t'. */
#define HAVE_PTRDIFF_T 1

/* Define to 1 if the system has the type `ssize_t'. */
#define HAVE_SSIZE_T 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if `st_blksize' is member of `struct stat'. */
#define HAVE_STRUCT_STAT_ST_BLKSIZE 1

/* Define to 1 if your `struct stat' has `st_blksize'. Deprecated, use
   `HAVE_STRUCT_STAT_ST_BLKSIZE' instead. */
#define HAVE_ST_BLKSIZE 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if the system has the type `uchar'. */
/* #undef HAVE_UCHAR */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Type of NC_BYTE */
#define NCBYTE_T byte

/* Type of NC_SHORT */
#define NCSHORT_T integer*2

/* C type for Fortran dobule */
/* #undef NF_DOUBLEPRECISION_IS_C_ */

/* C type for Fortran INT1 */
/* #undef NF_INT1_IS_C_ */

/* Type for Fortran INT1 */
#define NF_INT1_T byte

/* C type for Fortran INT2 */
/* #undef NF_INT2_IS_C_ */

/* Type for Fortran INT2 */
#define NF_INT2_T integer*2

/* C type for Fortran INT */
/* #undef NF_INT_IS_C_ */

/* C type for Fortran REAL */
/* #undef NF_REAL_IS_C_ */

/* Does sytem have IEEE FLOAT */
/* #undef NO_IEEE_FLOAT */

/* Define if system lacks strerror */
/* #undef NO_STRERROR */

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME ""

/* Define to the full name and version of this package. */
#define PACKAGE_STRING ""

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME ""

/* Define to the version of this package. */
#define PACKAGE_VERSION ""

/* full pnetcdf version string */
#define PNETCDF_VERSION "1.2.0"

/* major version number */
#define PNETCDF_VERSION_MAJOR 1

/* minor version number */
#define PNETCDF_VERSION_MINOR 2

/* pre-release string */
#define PNETCDF_VERSION_PRE 

/* sub version number */
#define PNETCDF_VERSION_SUB 0

/* The size of `double', as computed by sizeof. */
#define SIZEOF_DOUBLE 8

/* The size of `float', as computed by sizeof. */
#define SIZEOF_FLOAT 4

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* The number of bytes in an MPI_Offset */
#define SIZEOF_MPI_OFFSET 8

/* The size of `off_t', as computed by sizeof. */
#define SIZEOF_OFF_T 8

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT 2

/* The size of `size_t', as computed by sizeof. */
#define SIZEOF_SIZE_T 8

/* If using the C implementation of alloca, define if you know the
   direction of stack growth for your system; otherwise it will be
   automatically deduced at runtime.
	STACK_DIRECTION > 0 => grows toward higher addresses
	STACK_DIRECTION < 0 => grows toward lower addresses
	STACK_DIRECTION = 0 => direction of growth unknown */
/* #undef STACK_DIRECTION */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if your processor stores words with the most significant byte
   first (like Motorola and SPARC, unlike Intel and VAX). */
/* #undef WORDS_BIGENDIAN */

/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
#define YYTEXT_POINTER 1

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Define to 1 if type `char' is unsigned and you are not using gcc.  */
#ifndef __CHAR_UNSIGNED__
/* # undef __CHAR_UNSIGNED__ */
#endif

/* Define to `long int' if <sys/types.h> does not define. */
/* #undef off_t */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */
