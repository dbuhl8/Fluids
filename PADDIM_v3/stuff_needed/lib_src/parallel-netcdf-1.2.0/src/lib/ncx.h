/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* "$Id: ncx.h 891 2010-08-13 22:31:25Z robl $" */

#ifndef _NCX_H_
#define _NCX_H_

/*
 * An external data representation interface.
 *
 * This started out as a general replacement for ONC XDR,
 * specifically, the xdrmem family of functions.
 * 
 * We eventually realized that we could write more portable
 * code if we decoupled any association between the 'C' types
 * and the external types. (XDR has this association between the 'C'
 * types and the external representations, like xdr_int() takes
 * an int argument and goes to an external int representation.) 
 * So, now there is a matrix of functions.
 * 
 */

#include "ncconfig.h" /* output of 'configure' */
#include "rnd.h"
#include <stddef.h> 
#include <errno.h>
#include <sys/types.h> /* off_t */
#include <limits.h>
#include <float.h>
#include "nctypes.h"
#include <stdint.h>

#include <mpi.h>

#include "pnetcdf.h"


/*
 * External sizes of the primitive elements.
 */
#define X_SIZEOF_CHAR		1
#define X_SIZEOF_SHORT		2
#define X_SIZEOF_INT		4	/* xdr_int */
/*#if 0*/
#define X_SIZEOF_LONG		8	/* xdr_long_long */
/*#endif*/
#define X_SIZEOF_FLOAT		4
#define X_SIZEOF_DOUBLE		8

/*
 * For now, netcdf is limited to 32 bit sizes. 
 * If compiled with support for "large files" then 
 * netcdf will use a 64 bit off_t and it can then write a file 
 * using 64 bit offsets.
 *  see also X_SIZE_MAX, X_OFF_MAX below
 */
#define X_SIZEOF_OFF_T		(sizeof(off_t))
#define X_SIZEOF_SIZE_T		X_SIZEOF_INT

/*
 * limits of the external representation
 * we rely on ANSI-C defined constants in limits.h. Do any modern environments
 * not have these?
 */
#define X_SCHAR_MIN	SCHAR_MIN
#define X_SCHAR_MAX	SCHAR_MAX
#define X_UCHAR_MAX	UCHAR_MAX
#define X_SHORT_MIN	SHRT_MIN
#define X_SHRT_MIN	X_SHORT_MIN	/* alias compatible with limits.h */
#define X_SHORT_MAX	SHRT_MAX
#define X_SHRT_MAX	X_SHORT_MAX	/* alias compatible with limits.h */
#define X_USHORT_MAX	USHRT_MAX
#define X_USHRT_MAX	X_USHORT_MAX	/* alias compatible with limits.h */
#define X_INT_MIN	INT_MIN
#define X_INT_MAX	INT_MAX
#define X_UINT_MAX	UINT_MAX
#if 0
#define X_LONG_MIN	(-2147483647-1)
#define X_LONG_MAX	2147483647
#define X_ULONG_MAX	4294967295U
#endif
#define X_FLOAT_MAX	FLT_MAX
#define X_FLOAT_MIN	-(FLT_MAX)
#define X_FLT_MAX	X_FLOAT_MAX	/* alias compatible with limits.h */
/* scalb(1. - scalb(.5 , -52), 1024) */
#define X_DOUBLE_MAX	DBL_MAX
#define X_DOUBLE_MIN	-(DBL_MAX)
#define X_DBL_MAX	X_DOUBLE_MAX	/* alias compatible with limits.h */

#define X_SIZE_MAX	X_UINT_MAX	
#define X_OFF_MAX	X_INT_MAX

#define X_INT64_T_MAX	9223372036854775808
#define X_UINT64_T_MAX	18446744073709551616U

/* Begin ncx_len */

/*
 * ncx_len_xxx() interfaces are defined as macros below, 
 * These give the length of an array of nelems of the type.
 * N.B. The 'char' and 'short' interfaces give the X_ALIGNED length.
 */
#define X_ALIGN			4	/* a.k.a. BYTES_PER_XDR_UNIT */

#define ncmpix_len_char(nelems) \
	_RNDUP((nelems), X_ALIGN)

#define ncmpix_len_short(nelems) \
	(((nelems) + (nelems)%2)  * X_SIZEOF_SHORT)

#define ncmpix_len_int(nelems) \
	((nelems) * X_SIZEOF_INT)

#define ncmpix_len_long(nelems) \
	((nelems) * X_SIZEOF_LONG)

#define ncmpix_len_float(nelems) \
	((nelems) * X_SIZEOF_FLOAT)

#define ncmpix_len_double(nelems) \
	((nelems) * X_SIZEOF_DOUBLE)

/* End ncx_len */

#ifdef __CHAR_UNSIGNED__
	/* 'char' is unsigned, declare ncbyte as 'signed char' */
typedef signed char schar;

#else
	/* 'char' is signed */
typedef signed char schar;

#endif	/* __CHAR_UNSIGNED__ */

/*
 * Primitive numeric conversion functions.
 * The `put' functions convert from native internal
 * type to the external type, while the `get' functions
 * convert from the external to the internal.
 *
 * These take the form
 *	int ncx_get_{external_type}_{internal_type}(
 *		const void *xp,
 *		internal_type *ip
 *	);
 *	int ncx_put_{external_type}_{internal_type}(
 *		void *xp,
 *		const internal_type *ip
 *	);
 * where
 *	`external_type' and `internal_type' chosen from
		schar
		uchar
		short
		ushort
		int
		uint
		long
		ulong
		float
		double
 *
 * Not all combinations make sense.
 * We may not implement all combinations that make sense.
 * The netcdf functions that use this ncx interface don't
 * use these primitive conversion functions. They use the
 * aggregate conversion functions declared below.
 *
 * Storage for a single element of external type is at the `void * xp'
 * argument.
 *
 * Storage for a single element of internal type is at `ip' argument.
 *
 * These functions return 0 (ENOERR) when no error occured,
 * or NC_ERANGE when the value being converted is too large.
 * When NC_ERANGE occurs, an undefined (implementation dependent)
 * conversion may have occured.
 *
 * Note that loss of precision may occur silently.
 *
 */

#if 0
extern int
ncmpix_get_schar_schar(const void *xp, schar *ip);
extern int
ncmpix_get_schar_uchar(const void *xp, uchar *ip);
extern int
ncmpix_get_schar_short(const void *xp, short *ip);
extern int
ncmpix_get_schar_int(const void *xp, int *ip);
extern int
ncmpix_get_schar_long(const void *xp, long *ip);
extern int
ncmpix_get_schar_float(const void *xp, float *ip);
extern int
ncmpix_get_schar_double(const void *xp, double *ip);

extern int
ncmpix_put_schar_schar(void *xp, const schar *ip);
extern int
ncmpix_put_schar_uchar(void *xp, const uchar *ip);
extern int
ncmpix_put_schar_short(void *xp, const short *ip);
extern int
ncmpix_put_schar_int(void *xp, const int *ip);
extern int
ncmpix_put_schar_long(void *xp, const long *ip);
extern int
ncmpix_put_schar_float(void *xp, const float *ip);
extern int
ncmpix_put_schar_double(void *xp, const double *ip);
#endif
 

extern int
ncmpix_get_short_schar(const void *xp, schar *ip);
extern int
ncmpix_get_short_uchar(const void *xp, uchar *ip);
extern int
ncmpix_get_short_short(const void *xp, short *ip);
extern int
ncmpix_get_short_int(const void *xp, int *ip);
extern int
ncmpix_get_short_long(const void *xp, long *ip);
extern int
ncmpix_get_short_float(const void *xp, float *ip);
extern int
ncmpix_get_short_double(const void *xp, double *ip);

extern int
ncmpix_put_short_schar(void *xp, const schar *ip);
extern int
ncmpix_put_short_uchar(void *xp, const uchar *ip);
extern int
ncmpix_put_short_short(void *xp, const short *ip);
extern int
ncmpix_put_short_int(void *xp, const int *ip);
extern int
ncmpix_put_short_long(void *xp, const long *ip);
extern int
ncmpix_put_short_float(void *xp, const float *ip);
extern int
ncmpix_put_short_double(void *xp, const double *ip);
 

extern int
ncmpix_get_int_schar(const void *xp, schar *ip);
extern int
ncmpix_get_int_uchar(const void *xp, uchar *ip);
extern int
ncmpix_get_int_short(const void *xp, short *ip);
extern int
ncmpix_get_int_int(const void *xp, int *ip);
extern int
ncmpix_get_int_long(const void *xp, long *ip);
extern int
ncmpix_get_long_long(const void *xp,MPI_Offset *ip);
extern int
ncmpix_get_int_float(const void *xp, float *ip);
extern int
ncmpix_get_int_double(const void *xp, double *ip);

extern int
ncmpix_put_int_schar(void *xp, const schar *ip);
extern int
ncmpix_put_int_uchar(void *xp, const uchar *ip);
extern int
ncmpix_put_int_short(void *xp, const short *ip);
extern int
ncmpix_put_int_int(void *xp, const int *ip);
extern int
ncmpix_put_int_long(void *xp, const long *ip);
extern int
ncmpix_put_int_float(void *xp, const float *ip);
extern int
ncmpix_put_int_double(void *xp, const double *ip);
 

extern int
ncmpix_get_float_schar(const void *xp, schar *ip);
extern int
ncmpix_get_float_uchar(const void *xp, uchar *ip);
extern int
ncmpix_get_float_short(const void *xp, short *ip);
extern int
ncmpix_get_float_int(const void *xp, int *ip);
extern int
ncmpix_get_float_long(const void *xp, long *ip);
extern int
ncmpix_get_float_float(const void *xp, float *ip);
extern int
ncmpix_get_float_double(const void *xp, double *ip);

extern int
ncmpix_put_float_schar(void *xp, const schar *ip);
extern int
ncmpix_put_float_uchar(void *xp, const uchar *ip);
extern int
ncmpix_put_float_short(void *xp, const short *ip);
extern int
ncmpix_put_float_int(void *xp, const int *ip);
extern int
ncmpix_put_float_long(void *xp, const long *ip);
extern int
ncmpix_put_float_float(void *xp, const float *ip);
extern int
ncmpix_put_float_double(void *xp, const double *ip);
 

extern int
ncmpix_get_double_schar(const void *xp, schar *ip);
extern int
ncmpix_get_double_uchar(const void *xp, uchar *ip);
extern int
ncmpix_get_double_short(const void *xp, short *ip);
extern int
ncmpix_get_double_int(const void *xp, int *ip);
extern int
ncmpix_get_double_long(const void *xp, long *ip);
extern int
ncmpix_get_double_float(const void *xp, float *ip);
extern int
ncmpix_get_double_double(const void *xp, double *ip);

extern int
ncmpix_put_double_schar(void *xp, const schar *ip);
extern int
ncmpix_put_double_uchar(void *xp, const uchar *ip);
extern int
ncmpix_put_double_short(void *xp, const short *ip);
extern int
ncmpix_put_double_int(void *xp, const int *ip);
extern int
ncmpix_put_double_long(void *xp, const long *ip);
extern int
ncmpix_put_double_float(void *xp, const float *ip);
extern int
ncmpix_put_double_double(void *xp, const double *ip);
 

/*
 * Other primitive conversion functions
 * N.B. slightly different interface
 * Used by netcdf.
 */

/* ncx_get_int_size_t */
extern int
ncmpix_get_size_t(const void **xpp, MPI_Offset *ulp, MPI_Offset sizeof_off_t);
/* ncx_get_int_off_t */
extern int
ncmpix_get_off_t(const void **xpp, MPI_Offset *lp, MPI_Offset sizeof_off_t);

/* ncx_put_int_size_t */
extern int
ncmpix_put_size_t(void **xpp, const MPI_Offset *ulp, MPI_Offset sizeof_off_t);
extern int
ncmpix_put_size_t1(void **xpp, const MPI_Offset *ulp);
/* ncx_put_int_off_t */
extern int
ncmpix_put_off_t(void **xpp, const MPI_Offset *lp, MPI_Offset sizeof_off_t);


/*
 * Aggregate numeric conversion functions.
 * Convert an array.  Replaces xdr_array(...).
 * These functions are used by netcdf. Unlike the xdr
 * interface, we optimize for aggregate conversions.
 * This functions should be implemented to take advantage
 * of multiple processor / parallel hardware where available.
 *
 * These take the form
 *	int ncx_getn_{external_type}_{internal_type}(
 *		const void *xpp,
 *		size_t nelems,
 *		internal_type *ip
 *	);
 *	int ncx_putn_{external_type}_{internal_type}(
 *		void **xpp,
 *		size_t nelems,
 *		const internal_type *ip
 *	);
 * Where the types are as in the primitive numeric conversion functions.
 *
 * The value of the pointer to pointer argument, *xpp, is
 * expected to reference storage for `nelems' of the external
 * type.  On return, it modified to reference just past the last
 * converted external element.
 *
 * The types whose external size is less than X_ALIGN also have `pad'
 * interfaces. These round (and zero fill on put) *xpp up to X_ALIGN
 * boundaries. (This is the usual xdr behavior.)
 *
 * The `ip' argument should point to an array of `nelems' of
 * internal_type.
 *
 * Range errors (NC_ERANGE) for a individual values in the array 
 * DO NOT terminate the array conversion. All elements are converted,
 * with some having undefined values.
 * If any range error occurs, the function returns NC_ERANGE.
 *
 */

extern int
ncmpix_getn_schar_schar(const void **xpp, MPI_Offset nelems, schar *ip);
extern int
ncmpix_getn_schar_uchar(const void **xpp, MPI_Offset nelems, uchar *ip);
extern int
ncmpix_getn_schar_short(const void **xpp, MPI_Offset nelems, short *ip);
extern int
ncmpix_getn_schar_int(const void **xpp, MPI_Offset nelems, int *ip);
extern int
ncmpix_getn_schar_long(const void **xpp, MPI_Offset nelems, long *ip);
extern int
ncmpix_getn_schar_float(const void **xpp, MPI_Offset nelems, float *ip);
extern int
ncmpix_getn_schar_double(const void **xpp, MPI_Offset nelems, double *ip);

extern int
ncmpix_pad_getn_schar_schar(const void **xpp, MPI_Offset nelems, schar *ip);
extern int
ncmpix_pad_getn_schar_uchar(const void **xpp, MPI_Offset nelems, uchar *ip);
extern int
ncmpix_pad_getn_schar_short(const void **xpp, MPI_Offset nelems, short *ip);
extern int
ncmpix_pad_getn_schar_int(const void **xpp, MPI_Offset nelems, int *ip);
extern int
ncmpix_pad_getn_schar_long(const void **xpp, MPI_Offset nelems, long *ip);
extern int
ncmpix_pad_getn_schar_float(const void **xpp, MPI_Offset nelems, float *ip);
extern int
ncmpix_pad_getn_schar_double(const void **xpp, MPI_Offset nelems, double *ip);

extern int
ncmpix_putn_schar_schar(void **xpp, MPI_Offset nelems, const schar *ip);
extern int
ncmpix_putn_schar_uchar(void **xpp, MPI_Offset nelems, const uchar *ip);
extern int
ncmpix_putn_schar_short(void **xpp, MPI_Offset nelems, const short *ip);
extern int
ncmpix_putn_schar_int(void **xpp, MPI_Offset nelems, const int *ip);
extern int
ncmpix_putn_schar_long(void **xpp, MPI_Offset nelems, const long *ip);
extern int
ncmpix_putn_schar_float(void **xpp, MPI_Offset nelems, const float *ip);
extern int
ncmpix_putn_schar_double(void **xpp, MPI_Offset nelems, const double *ip);
 
extern int
ncmpix_pad_putn_schar_schar(void **xpp, MPI_Offset nelems, const schar *ip);
extern int
ncmpix_pad_putn_schar_uchar(void **xpp, MPI_Offset nelems, const uchar *ip);
extern int
ncmpix_pad_putn_schar_short(void **xpp, MPI_Offset nelems, const short *ip);
extern int
ncmpix_pad_putn_schar_int(void **xpp, MPI_Offset nelems, const int *ip);
extern int
ncmpix_pad_putn_schar_long(void **xpp, MPI_Offset nelems, const long *ip);
extern int
ncmpix_pad_putn_schar_float(void **xpp, MPI_Offset nelems, const float *ip);
extern int
ncmpix_pad_putn_schar_double(void **xpp, MPI_Offset nelems, const double *ip);


extern int
ncmpix_getn_short_schar(const void **xpp, MPI_Offset nelems, schar *ip);
extern int
ncmpix_getn_short_uchar(const void **xpp, MPI_Offset nelems, uchar *ip);
extern int
ncmpix_getn_short_short(const void **xpp, MPI_Offset nelems, short *ip);
extern int
ncmpix_getn_short_int(const void **xpp, MPI_Offset nelems, int *ip);
extern int
ncmpix_getn_short_long(const void **xpp, MPI_Offset nelems, long *ip);
extern int
ncmpix_getn_short_float(const void **xpp, MPI_Offset nelems, float *ip);
extern int
ncmpix_getn_short_double(const void **xpp, MPI_Offset nelems, double *ip);

extern int
ncmpix_pad_getn_short_schar(const void **xpp, MPI_Offset nelems, schar *ip);
extern int
ncmpix_pad_getn_short_uchar(const void **xpp, MPI_Offset nelems, uchar *ip);
extern int
ncmpix_pad_getn_short_short(const void **xpp, MPI_Offset nelems, short *ip);
extern int
ncmpix_pad_getn_short_int(const void **xpp, MPI_Offset nelems, int *ip);
extern int
ncmpix_pad_getn_short_long(const void **xpp, MPI_Offset nelems, long *ip);
extern int
ncmpix_pad_getn_short_float(const void **xpp, MPI_Offset nelems, float *ip);
extern int
ncmpix_pad_getn_short_double(const void **xpp, MPI_Offset nelems, double *ip);

extern int
ncmpix_putn_short_schar(void **xpp, MPI_Offset nelems, const schar *ip);
extern int
ncmpix_putn_short_uchar(void **xpp, MPI_Offset nelems, const uchar *ip);
extern int
ncmpix_putn_short_short(void **xpp, MPI_Offset nelems, const short *ip);
extern int
ncmpix_putn_short_int(void **xpp, MPI_Offset nelems, const int *ip);
extern int
ncmpix_putn_short_long(void **xpp, MPI_Offset nelems, const long *ip);
extern int
ncmpix_putn_short_float(void **xpp, MPI_Offset nelems, const float *ip);
extern int
ncmpix_putn_short_double(void **xpp, MPI_Offset nelems, const double *ip);
 
extern int
ncmpix_pad_putn_short_schar(void **xpp, MPI_Offset nelems, const schar *ip);
extern int
ncmpix_pad_putn_short_uchar(void **xpp, MPI_Offset nelems, const uchar *ip);
extern int
ncmpix_pad_putn_short_short(void **xpp, MPI_Offset nelems, const short *ip);
extern int
ncmpix_pad_putn_short_int(void **xpp, MPI_Offset nelems, const int *ip);
extern int
ncmpix_pad_putn_short_long(void **xpp, MPI_Offset nelems, const long *ip);
extern int
ncmpix_pad_putn_short_float(void **xpp, MPI_Offset nelems, const float *ip);
extern int
ncmpix_pad_putn_short_double(void **xpp, MPI_Offset nelems, const double *ip);


extern int
ncmpix_getn_int_schar(const void **xpp, MPI_Offset nelems, schar *ip);
extern int
ncmpix_getn_int_uchar(const void **xpp, MPI_Offset nelems, uchar *ip);
extern int
ncmpix_getn_int_short(const void **xpp, MPI_Offset nelems, short *ip);
extern int
ncmpix_getn_int_int(const void **xpp, MPI_Offset nelems, int *ip);
extern int
ncmpix_getn_int_long(const void **xpp, MPI_Offset nelems, long *ip);
extern int
ncmpix_getn_long_long(const void **xpp, MPI_Offset nelems, MPI_Offset *ip);
extern int
ncmpix_getn_int_float(const void **xpp, MPI_Offset nelems, float *ip);
extern int
ncmpix_getn_int_double(const void **xpp, MPI_Offset nelems, double *ip);

extern int
ncmpix_putn_int_schar(void **xpp, MPI_Offset nelems, const schar *ip);
extern int
ncmpix_putn_int_uchar(void **xpp, MPI_Offset nelems, const uchar *ip);
extern int
ncmpix_putn_int_short(void **xpp, MPI_Offset nelems, const short *ip);
extern int
ncmpix_putn_int_int(void **xpp, MPI_Offset nelems, const int *ip);
extern int
ncmpix_putn_int_long(void **xpp, MPI_Offset nelems, const long *ip);
extern int
ncmpix_putn_int_float(void **xpp, MPI_Offset nelems, const float *ip);
extern int
ncmpix_putn_int_double(void **xpp, MPI_Offset nelems, const double *ip);
 

extern int
ncmpix_getn_float_schar(const void **xpp, MPI_Offset nelems, schar *ip);
extern int
ncmpix_getn_float_uchar(const void **xpp, MPI_Offset nelems, uchar *ip);
extern int
ncmpix_getn_float_short(const void **xpp, MPI_Offset nelems, short *ip);
extern int
ncmpix_getn_float_int(const void **xpp, MPI_Offset nelems, int *ip);
extern int
ncmpix_getn_float_long(const void **xpp, MPI_Offset nelems, long *ip);
extern int
ncmpix_getn_float_float(const void **xpp, MPI_Offset nelems, float *ip);
extern int
ncmpix_getn_float_double(const void **xpp, MPI_Offset nelems, double *ip);

extern int
ncmpix_putn_float_schar(void **xpp, MPI_Offset nelems, const schar *ip);
extern int
ncmpix_putn_float_uchar(void **xpp, MPI_Offset nelems, const uchar *ip);
extern int
ncmpix_putn_float_short(void **xpp, MPI_Offset nelems, const short *ip);
extern int
ncmpix_putn_float_int(void **xpp, MPI_Offset nelems, const int *ip);
extern int
ncmpix_putn_float_long(void **xpp, MPI_Offset nelems, const long *ip);
extern int
ncmpix_putn_float_float(void **xpp, MPI_Offset nelems, const float *ip);
extern int
ncmpix_putn_float_double(void **xpp, MPI_Offset nelems, const double *ip);
 

extern int
ncmpix_getn_double_schar(const void **xpp, MPI_Offset nelems, schar *ip);
extern int
ncmpix_getn_double_uchar(const void **xpp, MPI_Offset nelems, uchar *ip);
extern int
ncmpix_getn_double_short(const void **xpp, MPI_Offset nelems, short *ip);
extern int
ncmpix_getn_double_int(const void **xpp, MPI_Offset nelems, int *ip);
extern int
ncmpix_getn_double_long(const void **xpp, MPI_Offset nelems, long *ip);
extern int
ncmpix_getn_double_float(const void **xpp, MPI_Offset nelems, float *ip);
extern int
ncmpix_getn_double_double(const void **xpp, MPI_Offset nelems, double *ip);

extern int
ncmpix_putn_double_schar(void **xpp, MPI_Offset nelems, const schar *ip);
extern int
ncmpix_putn_double_uchar(void **xpp, MPI_Offset nelems, const uchar *ip);
extern int
ncmpix_putn_double_short(void **xpp, MPI_Offset nelems, const short *ip);
extern int
ncmpix_putn_double_int(void **xpp, MPI_Offset nelems, const int *ip);
extern int
ncmpix_putn_double_long(void **xpp, MPI_Offset nelems, const long *ip);
extern int
ncmpix_putn_double_float(void **xpp, MPI_Offset nelems, const float *ip);
extern int
ncmpix_putn_double_double(void **xpp, MPI_Offset nelems, const double *ip);
 

/*
 * Other aggregate conversion functions.
 */

/* read ASCII characters */
extern int
ncmpix_getn_text(const void **xpp, MPI_Offset nchars, char *cp);
extern int
ncmpix_pad_getn_text(const void **xpp, MPI_Offset nchars, char *cp);

/* write ASCII characters */
extern int
ncmpix_putn_text(void **xpp, MPI_Offset nchars, const char *cp);
extern int
ncmpix_pad_putn_text(void **xpp, MPI_Offset nchars, const char *cp);

/* for symmetry */
#define ncmpix_getn_char_char(xpp, nelems, fillp) ncmpix_getn_text(xpp, nelems, fillp)
#define ncmpix_putn_char_char(xpp, nelems, fillp) ncmpix_putn_text(xpp, nelems, fillp)

/* read opaque data */
extern int
ncmpix_getn_void(const void **xpp, MPI_Offset nchars, void *vp);
extern int
ncmpix_pad_getn_void(const void **xpp, MPI_Offset nchars, void *vp);

/* write opaque data */
extern int
ncmpix_putn_void(void **xpp, MPI_Offset nchars, const void *vp);
extern int
ncmpix_pad_putn_void(void **xpp, MPI_Offset nchars, const void *vp);

#endif /* _NCX_H_ */
