/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *	Copyright 1996, University Corporation for Atmospheric Research
 *	See netcdf/COPYRIGHT file for copying and redistribution conditions.
 * 	
 * 	This file contains some routines derived from code
 *	which is copyrighted by Sun Microsystems, Inc.
 *	The "#ifdef vax" versions of
 *		 ncmpix_put_float_float()
 *		 ncmpix_get_float_float()
 *		 ncmpix_put_double_double()
 *		 ncmpix_get_double_double()
 *		 ncmpix_putn_float_float()
 *		 ncmpix_getn_float_float()
 *		 ncmpix_putn_double_double()
 *		 ncmpix_getn_double_double()
 * 	are derived from xdr_float() and xdr_double() routines
 *	in the freely available, copyrighted Sun RPCSRC 3.9
 *	distribution, xdr_float.c.
 * 	Our "value added" is that these are always memory to memory,
 *	they handle IEEE subnormals properly, and their "n" versions
 *	operate speedily on arrays.
 */
/* $Id: ncx.c 836 2010-06-01 22:22:05Z robl $ */

/*
 * An external data representation interface.
 */

#include "nc.h"
#include "ncx.h"
#include <string.h>
#include <limits.h>
#include <sys/types.h>
/* alias poorly named limits.h macros */
#define  SHORT_MAX  SHRT_MAX
#define  SHORT_MIN  SHRT_MIN
#define USHORT_MAX USHRT_MAX
#include <float.h>
#ifndef FLT_MAX /* This POSIX macro missing on some systems */
# ifndef NO_IEEE_FLOAT
# define FLT_MAX 3.40282347e+38f
# else
# error "You will need to define FLT_MAX"
# endif
#endif
#include <assert.h>

/*
 * If the machine's float domain is "smaller" than the external one
 * use the machine domain
 */
#if defined(FLT_MAX_EXP) && FLT_MAX_EXP < 128 /* 128 is X_FLT_MAX_EXP */
#undef X_FLOAT_MAX
# define X_FLOAT_MAX FLT_MAX
#undef X_FLOAT_MIN
# define X_FLOAT_MIN (-X_FLOAT_MAX)
#endif

#if defined(_SX) && _SX != 0 /* NEC SUPER UX */
#if _INT64
#undef  INT_MAX /* workaround cpp bug */
#define INT_MAX  X_INT_MAX
#undef  INT_MIN /* workaround cpp bug */
#define INT_MIN  X_INT_MIN
#undef  LONG_MAX /* workaround cpp bug */
#define LONG_MAX  X_INT_MAX
#undef  LONG_MIN /* workaround cpp bug */
#define LONG_MIN  X_INT_MIN
#elif _LONG64
#undef  LONG_MAX /* workaround cpp bug */
#define LONG_MAX  4294967295L
#undef  LONG_MIN /* workaround cpp bug */
#define LONG_MIN -4294967295L
#endif
#endif /* _SX */

static const char nada[X_ALIGN] = {0, 0, 0, 0};

#ifndef WORDS_BIGENDIAN
/* LITTLE_ENDIAN: DEC and intel */
/*
 * Routines to convert to BIGENDIAN.
 * Optimize the swapn?b() and swap?b() routines aggressivly.
 */

#define SWAP2(a) ( (((a) & 0xff) << 8) | \
		(((a) >> 8) & 0xff) )

#define SWAP4(a) ( ((a) << 24) | \
		(((a) <<  8) & 0x00ff0000) | \
		(((a) >>  8) & 0x0000ff00) | \
		(((a) >> 24) & 0x000000ff) )

/* netcdf-3.6.2beta5 added loop unrolling to many of these routines.  Could
 * confer a 22% performance increase on little endian platforms if compiler
 * does not already aggressively unroll loops */

static void
swapn2b(void *dst, const void *src, MPI_Offset nn)
{
	char *op = dst;
	const char *ip = src;

	while(nn > 3)
	{
		*op++ = *(++ip);
		*op++ = *(ip++ -1);
		*op++ = *(++ip);
		*op++ = *(ip++ -1);
		*op++ = *(++ip);
		*op++ = *(ip++ -1);
		*op++ = *(++ip);
		*op++ = *(ip++ -1);
		nn -= 4;
	}
	while(nn-- != 0)
	{
		*op++ = *(++ip);
		*op++ = *(ip++ -1);
	}
}

# ifndef vax
static void
swap4b(void *dst, const void *src)
{
	char *op = dst;
	const char *ip = src;
	op[0] = ip[3];
	op[1] = ip[2];
	op[2] = ip[1];
	op[3] = ip[0];
}
# endif /* !vax */

static void
swapn4b(void *dst, const void *src, MPI_Offset nn)
{
	char *op = dst;
	const char *ip = src;

	while(nn > 3)
	{
		op[0] = ip[3];
		op[1] = ip[2];
		op[2] = ip[1];
		op[3] = ip[0];
		op[4] = ip[7];
		op[5] = ip[6];
		op[6] = ip[5];
		op[7] = ip[4];
		op[8] = ip[11];
		op[9] = ip[10];
		op[10] = ip[9];
		op[11] = ip[8];
		op[12] = ip[15];
		op[13] = ip[14];
		op[14] = ip[13];
		op[15] = ip[12];
		op += 16;
		ip += 16;
		nn -= 4;
	}
	while(nn-- != 0)
	{
		op[0] = ip[3];
		op[1] = ip[2];
		op[2] = ip[1];
		op[3] = ip[0];
		op += 4;
		ip += 4;
	}
}

# ifndef vax
static void
swap8b(void *dst, const void *src)
{
	char *op = dst;
	const char *ip = src;
	op[0] = ip[7];
	op[1] = ip[6];
	op[2] = ip[5];
	op[3] = ip[4];
	op[4] = ip[3];
	op[5] = ip[2];
	op[6] = ip[1];
	op[7] = ip[0];
}
# endif /* !vax */

# ifndef vax
static void
swapn8b(void *dst, const void *src, MPI_Offset nn)
{
	char *op = dst;
	const char *ip = src;

	while(nn > 1)
	{
		op[0] = ip[7];
		op[1] = ip[6];
		op[2] = ip[5];
		op[3] = ip[4];
		op[4] = ip[3];
		op[5] = ip[2];
		op[6] = ip[1];
		op[7] = ip[0];
		op[8] = ip[15];
		op[9] = ip[14];
		op[10] = ip[13];
		op[11] = ip[12];
		op[12] = ip[11];
		op[13] = ip[10];
		op[14] = ip[9];
		op[15] = ip[8];
		op += 16;
		ip += 16;
		nn -= 2;
	}
	while(nn-- != 0)
	{
		op[0] = ip[7];
		op[1] = ip[6];
		op[2] = ip[5];
		op[3] = ip[4];
		op[4] = ip[3];
		op[5] = ip[2];
		op[6] = ip[1];
		op[7] = ip[0];
		op += 8;
		ip += 8;
	}
}
# endif /* !vax */

#endif /* LITTLE_ENDIAN */


/*
 * Primitive numeric conversion functions.
 */

/* x_schar */

 /* We don't implement and x_schar primitives. */


/* x_short */

#if SHORT_MAX == X_SHORT_MAX
typedef short ix_short;
#define SIZEOF_IX_SHORT SIZEOF_SHORT
#define IX_SHORT_MAX SHORT_MAX
#elif INT_MAX >= X_SHORT_MAX
typedef int ix_short;
#define SIZEOF_IX_SHORT SIZEOF_INT
#define IX_SHORT_MAX INT_MAX
#elif LONG_MAX >= X_SHORT_MAX
typedef long ix_short;
#define SIZEOF_IX_SHORT SIZEOF_LONG
#define IX_SHORT_MAX LONG_MAX
#else
#error "ix_short implementation"
#endif

static void
get_ix_short(const void *xp, ix_short *ip)
{
	const uchar *cp = (const uchar *) xp;
	*ip = *cp++ << 8;
#if SIZEOF_IX_SHORT > X_SIZEOF_SHORT
	if(*ip & 0x8000)
	{
		/* extern is negative */
		*ip |= (~(0xffff)); /* N.B. Assumes "twos complement" */
	}
#endif
	*ip |= *cp; 
}

static void
put_ix_short(void *xp, const ix_short *ip)
{
	uchar *cp = (uchar *) xp;
	*cp++ = (*ip) >> 8;
	*cp = (*ip) & 0xff;
}


int
ncmpix_get_short_schar(const void *xp, schar *ip)
{
	ix_short xx=0;
	get_ix_short(xp, &xx);
	*ip = xx;
	if(xx > SCHAR_MAX || xx < SCHAR_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_short_uchar(const void *xp, uchar *ip)
{
	ix_short xx=0;
	get_ix_short(xp, &xx);
	*ip = xx;
	if(xx > UCHAR_MAX || xx < 0)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_short_short(const void *xp, short *ip)
{
#if SIZEOF_IX_SHORT == SIZEOF_SHORT && IX_SHORT_MAX == SHORT_MAX
	get_ix_short(xp, (ix_short *)ip);
	return NC_NOERR;
#else
	ix_short xx;
	get_ix_short(xp, &xx);
	*ip = xx;
#   if IX_SHORT_MAX > SHORT_MAX
	if(xx > SHORT_MAX || xx < SHORT_MIN)
		return NC_ERANGE;
#   endif
	return NC_NOERR;
#endif
}

int
ncmpix_get_short_int(const void *xp, int *ip)
{
#if SIZEOF_IX_SHORT == SIZEOF_INT && IX_SHORT_MAX == INT_MAX
	get_ix_short(xp, (ix_short *)ip);
	return NC_NOERR;
#else
	ix_short xx;
	get_ix_short(xp, &xx);
	*ip = xx;
#   if IX_SHORT_MAX > INT_MAX
	if(xx > INT_MAX || xx < INT_MIN)
		return NC_ERANGE;
#   endif
	return NC_NOERR;
#endif
}

int
ncmpix_get_short_long(const void *xp, long *ip)
{
#if SIZEOF_IX_SHORT == SIZEOF_LONG && IX_SHORT_MAX == LONG_MAX
	get_ix_short(xp, (ix_short *)ip);
	return NC_NOERR;
#else
	/* assert(LONG_MAX >= X_SHORT_MAX); */
	ix_short xx;
	get_ix_short(xp, &xx);
	*ip = xx;
	return NC_NOERR;
#endif
}

int
ncmpix_get_short_float(const void *xp, float *ip)
{
	ix_short xx;
	get_ix_short(xp, &xx);
	*ip = xx;
#if 0	/* TODO: determine when necessary */
	if(xx > FLT_MAX || xx < (-FLT_MAX))
		return NC_ERANGE;
#endif
	return NC_NOERR;
}

int
ncmpix_get_short_double(const void *xp, double *ip)
{
	/* assert(DBL_MAX >= X_SHORT_MAX); */
	ix_short xx;
	get_ix_short(xp, &xx);
	*ip = xx;
	return NC_NOERR;
}

int
ncmpix_put_short_schar(void *xp, const schar *ip)
{
	uchar *cp = (uchar *) xp;
	if(*ip & 0x80)
		*cp++ = 0xff;
	else
		*cp++ = 0;
	*cp = (uchar)*ip;
	return NC_NOERR;
}

int
ncmpix_put_short_uchar(void *xp, const uchar *ip)
{
	uchar *cp = (uchar *) xp;
	*cp++ = 0;
	*cp = *ip;
	return NC_NOERR;
}

int
ncmpix_put_short_short(void *xp, const short *ip)
{
#if SIZEOF_IX_SHORT == SIZEOF_SHORT && X_SHORT_MAX == SHORT_MAX
	put_ix_short(xp, (const ix_short *)ip);
	return NC_NOERR;
#else
	ix_short xx = (ix_short)*ip;
	put_ix_short(xp, &xx);
# if X_SHORT_MAX < SHORT_MAX
	if(*ip > X_SHORT_MAX || *ip < X_SHORT_MIN)
		return NC_ERANGE;
# endif
	return NC_NOERR;
#endif
}

int
ncmpix_put_short_int(void *xp, const int *ip)
{
#if SIZEOF_IX_SHORT == SIZEOF_INT && X_SHORT_MAX == INT_MAX
	put_ix_short(xp, (const ix_short *)ip);
	return NC_NOERR;
#else
	ix_short xx = (ix_short)*ip;
	put_ix_short(xp, &xx);
# if X_SHORT_MAX < INT_MAX
	if(*ip > X_SHORT_MAX || *ip < X_SHORT_MIN)
		return NC_ERANGE;
# endif
	return NC_NOERR;
#endif
}

int
ncmpix_put_short_long(void *xp, const long *ip)
{
#if SIZEOF_IX_SHORT == SIZEOF_LONG && X_SHORT_MAX == LONG_MAX
	put_ix_short(xp, (const ix_short *)ip);
	return NC_NOERR;
#else
	ix_short xx = (ix_short)*ip;
	put_ix_short(xp, &xx);
# if X_SHORT_MAX < LONG_MAX
	if(*ip > X_SHORT_MAX || *ip < X_SHORT_MIN)
		return NC_ERANGE;
# endif
	return NC_NOERR;
#endif
}

int
ncmpix_put_short_float(void *xp, const float *ip)
{
	ix_short xx = *ip;
	put_ix_short(xp, &xx);
	if(*ip > X_SHORT_MAX || *ip < X_SHORT_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_put_short_double(void *xp, const double *ip)
{
	ix_short xx = *ip;
	put_ix_short(xp, &xx);
	if(*ip > X_SHORT_MAX || *ip < X_SHORT_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

/* x_int */

#if SHORT_MAX == X_INT_MAX
typedef short ix_int;
#define SIZEOF_IX_INT SIZEOF_SHORT
#define IX_INT_MAX SHORT_MAX
#elif INT_MAX  >= X_INT_MAX
typedef int ix_int;
#define SIZEOF_IX_INT SIZEOF_INT
#define IX_INT_MAX INT_MAX
#elif LONG_MAX  >= X_INT_MAX
typedef long ix_int;
#define SIZEOF_IX_INT SIZEOF_LONG
#define IX_INT_MAX LONG_MAX
#else
#error "ix_int implementation"
#endif


static void
get_ix_int(const void *xp, ix_int *ip)
{
	const uchar *cp = (const uchar *) xp;

	*ip = *cp++ << 24;
#if SIZEOF_IX_INT > X_SIZEOF_INT
	if(*ip & 0x80000000)
	{
		/* extern is negative */
		*ip |= (~(0xffffffff)); /* N.B. Assumes "twos complement" */
	}
#endif
	*ip |= (*cp++ << 16);
	*ip |= (*cp++ << 8);
	*ip |= *cp; 
}

static void
put_ix_int(void *xp, const ix_int *ip)
{
	uchar *cp = (uchar *) xp;

	*cp++ = (*ip) >> 24;
	*cp++ = ((*ip) & 0x00ff0000) >> 16;
	*cp++ = ((*ip) & 0x0000ff00) >>  8;
	*cp   = ((*ip) & 0x000000ff);
}


int
ncmpix_get_int_schar(const void *xp, schar *ip)
{
	ix_int xx=0;
	get_ix_int(xp, &xx);
	*ip = xx;
	if(xx > SCHAR_MAX || xx < SCHAR_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_int_uchar(const void *xp, uchar *ip)
{
	ix_int xx=0;
	get_ix_int(xp, &xx);
	*ip = xx;
	if(xx > UCHAR_MAX || xx < 0)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_int_short(const void *xp, short *ip)
{
#if SIZEOF_IX_INT == SIZEOF_SHORT && IX_INT_MAX == SHORT_MAX
	get_ix_int(xp, (ix_int *)ip);
	return NC_NOERR;
#else
	ix_int xx=0;
	get_ix_int(xp, &xx);
	*ip = xx;
#  if IX_INT_MAX > SHORT_MAX
	if(xx > SHORT_MAX || xx < SHORT_MIN)
		return NC_ERANGE;
#  endif
	return NC_NOERR;
#endif
}

int
ncmpix_get_int_int(const void *xp, int *ip)
{
#if SIZEOF_IX_INT == SIZEOF_INT && IX_INT_MAX == INT_MAX
	get_ix_int(xp, (ix_int *)ip);
	return NC_NOERR;
#else
	ix_int xx;
	get_ix_int(xp, &xx);
	*ip = xx;
#  if IX_INT_MAX > INT_MAX
	if(xx > INT_MAX || xx < INT_MIN)
		return NC_ERANGE;
#  endif
	return NC_NOERR;
#endif
}

int
ncmpix_get_int_long(const void *xp, long *ip)
{
#if SIZEOF_IX_INT == SIZEOF_LONG && IX_INT_MAX == LONG_MAX
	get_ix_int(xp, (ix_int *)ip);
	return NC_NOERR;
#else
	ix_int xx;
	get_ix_int(xp, &xx);
	*ip = xx;
#  if IX_INT_MAX > LONG_MAX	/* unlikely */
	if(xx > LONG_MAX || xx < LONG_MIN)
		return NC_ERANGE;
#  endif
	return NC_NOERR;
#endif
}

int
ncmpix_get_long_long(const void *xp,MPI_Offset *ip)
{
#if SIZEOF_IX_INT == SIZEOF_LONG && IX_INT_MAX == LONG_MAX
	get_ix_int(xp, (ix_int *)ip);
	return NC_NOERR;
#else
	ix_int xx;
	get_ix_int(xp, &xx);
	*ip = xx;
#  if IX_INT_MAX > LONG_MAX	/* unlikely */
	if(xx > LONG_MAX || xx < LONG_MIN)
		return NC_ERANGE;
#  endif
	return NC_NOERR;
#endif
}

int
ncmpix_get_int_float(const void *xp, float *ip)
{
	ix_int xx;
	get_ix_int(xp, &xx);
	*ip = xx;
#if 0	/* TODO: determine when necessary */
	if(xx > FLT_MAX || xx < (-FLT_MAX))
		return NC_ERANGE;
#endif
	return NC_NOERR;
}

int
ncmpix_get_int_double(const void *xp, double *ip)
{
	/* assert((DBL_MAX >= X_INT_MAX); */
	ix_int xx;
	get_ix_int(xp, &xx);
	*ip = xx;
	return NC_NOERR;
}

int
ncmpix_put_int_schar(void *xp, const schar *ip)
{
	uchar *cp = (uchar *) xp;
	if(*ip & 0x80)
	{
		*cp++ = 0xff;
		*cp++ = 0xff;
		*cp++ = 0xff;
	}
	else
	{
		*cp++ = 0x00;
		*cp++ = 0x00;
		*cp++ = 0x00;
	}
	*cp = (uchar)*ip;
	return NC_NOERR;
}

int
ncmpix_put_int_uchar(void *xp, const uchar *ip)
{
	uchar *cp = (uchar *) xp;
	*cp++ = 0x00;
	*cp++ = 0x00;
	*cp++ = 0x00;
	*cp   = *ip;
	return NC_NOERR;
}

int
ncmpix_put_int_short(void *xp, const short *ip)
{
#if SIZEOF_IX_INT == SIZEOF_SHORT && IX_INT_MAX == SHORT_MAX
	put_ix_int(xp, (ix_int *)ip);
	return NC_NOERR;
#else
	ix_int xx = (ix_int)(*ip);
	put_ix_int(xp, &xx);
#   if IX_INT_MAX < SHORT_MAX
	if(*ip > X_INT_MAX || *ip < X_INT_MIN)
		return NC_ERANGE;
#   endif
	return NC_NOERR;
#endif
}

int
ncmpix_put_int_int(void *xp, const int *ip)
{
#if SIZEOF_IX_INT == SIZEOF_INT && IX_INT_MAX == INT_MAX
	put_ix_int(xp, (ix_int *)ip);
	return NC_NOERR;
#else
	ix_int xx = (ix_int)(*ip);
	put_ix_int(xp, &xx);
#   if IX_INT_MAX < INT_MAX
	if(*ip > X_INT_MAX || *ip < X_INT_MIN)
		return NC_ERANGE;
#   endif
	return NC_NOERR;
#endif
}

int
ncmpix_put_int_long(void *xp, const long *ip)
{
#if SIZEOF_IX_INT == SIZEOF_LONG && IX_INT_MAX == LONG_MAX
	put_ix_int(xp, (ix_int *)ip);
	return NC_NOERR;
#else
	ix_int xx = (ix_int)(*ip);
	put_ix_int(xp, &xx);
#   if IX_INT_MAX < LONG_MAX
	if(*ip > X_INT_MAX || *ip < X_INT_MIN)
		return NC_ERANGE;
#   endif
	return NC_NOERR;
#endif
}

int
ncmpix_put_int_float(void *xp, const float *ip)
{
	ix_int xx = (ix_int)(*ip);
	put_ix_int(xp, &xx);
	if(*ip > (double)X_INT_MAX || *ip < (double)X_INT_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_put_int_double(void *xp, const double *ip)
{
	ix_int xx = (ix_int)(*ip);
	put_ix_int(xp, &xx);
	if(*ip > X_INT_MAX || *ip < X_INT_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}
 

/* x_float */

#if X_SIZEOF_FLOAT == SIZEOF_FLOAT && !defined(NO_IEEE_FLOAT)

static void
get_ix_float(const void *xp, float *ip)
{
#ifdef WORDS_BIGENDIAN
	(void) memcpy(ip, xp, sizeof(float));
#else
	swap4b(ip, xp);
#endif
}

static void
put_ix_float(void *xp, const float *ip)
{
#ifdef WORDS_BIGENDIAN
	(void) memcpy(xp, ip, X_SIZEOF_FLOAT);
#else
	swap4b(xp, ip);
#endif
}

#elif vax

/* What IEEE single precision floating point looks like on a Vax */
struct	ieee_single {
	unsigned int	exp_hi       : 7;
	unsigned int	sign         : 1;
	unsigned int 	mant_hi      : 7;
	unsigned int	exp_lo       : 1;
	unsigned int	mant_lo_hi   : 8;
	unsigned int	mant_lo_lo   : 8;
};

/* Vax single precision floating point */
struct	vax_single {
	unsigned int	mantissa1 : 7;
	unsigned int	exp       : 8;
	unsigned int	sign      : 1;
	unsigned int	mantissa2 : 16;
};

#define VAX_SNG_BIAS	0x81
#define IEEE_SNG_BIAS	0x7f

static struct sgl_limits {
	struct vax_single s;
	struct ieee_single ieee;
} max = {
	{ 0x7f, 0xff, 0x0, 0xffff },	/* Max Vax */
	{ 0x7f, 0x0, 0x0, 0x1, 0x0, 0x0 }		/* Max IEEE */
};
static struct sgl_limits min = {
	{ 0x0, 0x0, 0x0, 0x0 },	/* Min Vax */
	{ 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 }		/* Min IEEE */
};

static void
get_ix_float(const void *xp, float *ip)
{
		struct vax_single *const vsp = (struct vax_single *) ip;
		const struct ieee_single *const isp =
			 (const struct ieee_single *) xp;
		unsigned exp = isp->exp_hi << 1 | isp->exp_lo;

		switch(exp) {
		case 0 :
			/* ieee subnormal */
			if(isp->mant_hi == min.ieee.mant_hi
				&& isp->mant_lo_hi == min.ieee.mant_lo_hi
				&& isp->mant_lo_lo == min.ieee.mant_lo_lo)
			{
				*vsp = min.s;
			}
			else
			{
				unsigned mantissa = (isp->mant_hi << 16)
					 | isp->mant_lo_hi << 8
					 | isp->mant_lo_lo;
				unsigned tmp = mantissa >> 20;
				if(tmp >= 4) {
					vsp->exp = 2;
				} else if (tmp >= 2) {
					vsp->exp = 1;
				} else {
					*vsp = min.s;
					break;
				} /* else */
				tmp = mantissa - (1 << (20 + vsp->exp ));
				tmp <<= 3 - vsp->exp;
				vsp->mantissa2 = tmp;
				vsp->mantissa1 = (tmp >> 16);
			}
			break;
		case 0xfe :
		case 0xff :
			*vsp = max.s;
			break;
		default :
			vsp->exp = exp - IEEE_SNG_BIAS + VAX_SNG_BIAS;
			vsp->mantissa2 = isp->mant_lo_hi << 8 | isp->mant_lo_lo;
			vsp->mantissa1 = isp->mant_hi;
		}

		vsp->sign = isp->sign;

}


static void
put_ix_float(void *xp, const float *ip)
{
		const struct vax_single *const vsp =
			 (const struct vax_single *)ip;
		struct ieee_single *const isp = (struct ieee_single *) xp;

		switch(vsp->exp){
		case 0 :
			/* all vax float with zero exponent map to zero */
			*isp = min.ieee;
			break;
		case 2 :
		case 1 :
		{
			/* These will map to subnormals */
			unsigned mantissa = (vsp->mantissa1 << 16)
					 | vsp->mantissa2;
			mantissa >>= 3 - vsp->exp;
			mantissa += (1 << (20 + vsp->exp));
			isp->mant_lo_lo = mantissa;
			isp->mant_lo_hi = mantissa >> 8;
			isp->mant_hi = mantissa >> 16;
			isp->exp_lo = 0;
			isp->exp_hi = 0;
		}
			break;
		case 0xff : /* max.s.exp */
			if( vsp->mantissa2 == max.s.mantissa2
				&& vsp->mantissa1 == max.s.mantissa1)
			{
				/* map largest vax float to ieee infinity */
				*isp = max.ieee;
				break;
			} /* else, fall thru */
		default :
		{
			unsigned exp = vsp->exp - VAX_SNG_BIAS + IEEE_SNG_BIAS;
			isp->exp_hi = exp >> 1;
			isp->exp_lo = exp;
			isp->mant_lo_lo = vsp->mantissa2;
			isp->mant_lo_hi = vsp->mantissa2 >> 8;
			isp->mant_hi = vsp->mantissa1;
		}
		}

		isp->sign = vsp->sign;

}

	/* vax */
#else
#error "ix_float implementation"
#endif


int
ncmpix_get_float_schar(const void *xp, schar *ip)
{
	float xx=0.0;
	get_ix_float(xp, &xx);
	*ip = (schar) xx;
	if(xx > SCHAR_MAX || xx < SCHAR_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_float_uchar(const void *xp, uchar *ip)
{
	float xx=0.0;
	get_ix_float(xp, &xx);
	*ip = (uchar) xx;
	if(xx > UCHAR_MAX || xx < 0)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_float_short(const void *xp, short *ip)
{
	float xx=0.0;
	get_ix_float(xp, &xx);
	*ip = (short) xx;
	if(xx > SHORT_MAX || xx < SHORT_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_float_int(const void *xp, int *ip)
{
	float xx=0.0;
	get_ix_float(xp, &xx);
	*ip = (int) xx;
	if(xx > (double)INT_MAX || xx < (double)INT_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_float_long(const void *xp, long *ip)
{
	float xx;
	get_ix_float(xp, &xx);
	*ip = (long) xx;
	if(xx > LONG_MAX || xx < LONG_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_float_float(const void *xp, float *ip)
{
	/* TODO */
	get_ix_float(xp, ip);
	return NC_NOERR;
}

int
ncmpix_get_float_double(const void *xp, double *ip)
{
	/* TODO */
	float xx;
	get_ix_float(xp, &xx);
	*ip = xx;
	return NC_NOERR;
}


int
ncmpix_put_float_schar(void *xp, const schar *ip)
{
	float xx = (float) *ip;
	put_ix_float(xp, &xx);
	return NC_NOERR;
}

int
ncmpix_put_float_uchar(void *xp, const uchar *ip)
{
	float xx = (float) *ip;
	put_ix_float(xp, &xx);
	return NC_NOERR;
}

int
ncmpix_put_float_short(void *xp, const short *ip)
{
	float xx = (float) *ip;
	put_ix_float(xp, &xx);
#if 0	/* TODO: figure this out */
	if((float)(*ip) > X_FLOAT_MAX || (float)(*ip) < X_FLOAT_MIN)
		return NC_ERANGE;
#endif
	return NC_NOERR;
}

int
ncmpix_put_float_int(void *xp, const int *ip)
{
	float xx = (float) *ip;
	put_ix_float(xp, &xx);
#if 1	/* TODO: figure this out */
	if((float)(*ip) > X_FLOAT_MAX || (float)(*ip) < X_FLOAT_MIN)
		return NC_ERANGE;
#endif
	return NC_NOERR;
}

int
ncmpix_put_float_long(void *xp, const long *ip)
{
	float xx = (float) *ip;
	put_ix_float(xp, &xx);
#if 1	/* TODO: figure this out */
	if((float)(*ip) > X_FLOAT_MAX || (float)(*ip) < X_FLOAT_MIN)
		return NC_ERANGE;
#endif
	return NC_NOERR;
}

int
ncmpix_put_float_float(void *xp, const float *ip)
{
	put_ix_float(xp, ip);
#ifdef NO_IEEE_FLOAT
	if(*ip > X_FLOAT_MAX || *ip < X_FLOAT_MIN)
		return NC_ERANGE;
#endif
	return NC_NOERR;
}

int
ncmpix_put_float_double(void *xp, const double *ip)
{
	float xx = (float) *ip;
	put_ix_float(xp, &xx);
	if(*ip > X_FLOAT_MAX || *ip < X_FLOAT_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

/* x_double */

#if X_SIZEOF_DOUBLE == SIZEOF_DOUBLE  && !defined(NO_IEEE_FLOAT)

static void
get_ix_double(const void *xp, double *ip)
{
#ifdef WORDS_BIGENDIAN
	(void) memcpy(ip, xp, sizeof(double));
#else
	swap8b(ip, xp);
#endif
}

static void
put_ix_double(void *xp, const double *ip)
{
#ifdef WORDS_BIGENDIAN
	(void) memcpy(xp, ip, X_SIZEOF_DOUBLE);
#else
	swap8b(xp, ip);
#endif
}

#elif vax

/* What IEEE double precision floating point looks like on a Vax */
struct	ieee_double {
	unsigned int	exp_hi   : 7;
	unsigned int	sign     : 1;
	unsigned int 	mant_6   : 4;
	unsigned int	exp_lo   : 4;
	unsigned int	mant_5   : 8;
	unsigned int	mant_4   : 8;

	unsigned int	mant_lo  : 32;
};

/* Vax double precision floating point */
struct  vax_double {
	unsigned int	mantissa1 : 7;
	unsigned int	exp       : 8;
	unsigned int	sign      : 1;
	unsigned int	mantissa2 : 16;
	unsigned int	mantissa3 : 16;
	unsigned int	mantissa4 : 16;
};

#define VAX_DBL_BIAS	0x81
#define IEEE_DBL_BIAS	0x3ff
#define MASK(nbits)	((1 << nbits) - 1)

static const struct dbl_limits {
	struct	vax_double d;
	struct	ieee_double ieee;
} dbl_limits[2] = {
	{{ 0x7f, 0xff, 0x0, 0xffff, 0xffff, 0xffff },	/* Max Vax */
	{ 0x7f, 0x0, 0x0, 0xf, 0x0, 0x0, 0x0}}, /* Max IEEE */
	{{ 0x0, 0x0, 0x0, 0x0, 0x0, 0x0},		/* Min Vax */
	{ 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0}}, /* Min IEEE */
};


static void
get_ix_double(const void *xp, double *ip)
{
	struct vax_double *const vdp =
			 (struct vax_double *)ip;
	const struct ieee_double *const idp =
			 (const struct ieee_double *) xp;
	{
		const struct dbl_limits *lim;
		int ii;
		for (ii = 0, lim = dbl_limits;
			ii < sizeof(dbl_limits)/sizeof(struct dbl_limits);
			ii++, lim++)
		{
			if ((idp->mant_lo == lim->ieee.mant_lo)
				&& (idp->mant_4 == lim->ieee.mant_4)
				&& (idp->mant_5 == lim->ieee.mant_5)
				&& (idp->mant_6 == lim->ieee.mant_6)
				&& (idp->exp_lo == lim->ieee.exp_lo)
				&& (idp->exp_hi == lim->ieee.exp_hi)
				)
			{
				*vdp = lim->d;
				goto doneit;
			}
		}
	}
	{
		unsigned exp = idp->exp_hi << 4 | idp->exp_lo;
		vdp->exp = exp - IEEE_DBL_BIAS + VAX_DBL_BIAS;
	}
	{
		unsigned mant_hi = ((idp->mant_6 << 16)
				 | (idp->mant_5 << 8)
				 | idp->mant_4);
		unsigned mant_lo = SWAP4(idp->mant_lo);
		vdp->mantissa1 = (mant_hi >> 13);
		vdp->mantissa2 = ((mant_hi & MASK(13)) << 3)
				| (mant_lo >> 29);
		vdp->mantissa3 = (mant_lo >> 13);
		vdp->mantissa4 = (mant_lo << 3);
	}
	doneit:
		vdp->sign = idp->sign;

}


static void
put_ix_double(void *xp, const double *ip)
{
	const struct vax_double *const vdp = 
			(const struct vax_double *)ip;
	struct ieee_double *const idp =
			 (struct ieee_double *) xp;

	if ((vdp->mantissa4 > (dbl_limits[0].d.mantissa4 - 3)) &&
		(vdp->mantissa3 == dbl_limits[0].d.mantissa3) &&
		(vdp->mantissa2 == dbl_limits[0].d.mantissa2) &&
		(vdp->mantissa1 == dbl_limits[0].d.mantissa1) &&
		(vdp->exp == dbl_limits[0].d.exp))
	{
		*idp = dbl_limits[0].ieee;
		goto shipit;
	}
	if ((vdp->mantissa4 == dbl_limits[1].d.mantissa4) &&
		(vdp->mantissa3 == dbl_limits[1].d.mantissa3) &&
		(vdp->mantissa2 == dbl_limits[1].d.mantissa2) &&
		(vdp->mantissa1 == dbl_limits[1].d.mantissa1) &&
		(vdp->exp == dbl_limits[1].d.exp))
	{
		*idp = dbl_limits[1].ieee;
		goto shipit;
	}

	{
		unsigned exp = vdp->exp - VAX_DBL_BIAS + IEEE_DBL_BIAS;

		unsigned mant_lo = ((vdp->mantissa2 & MASK(3)) << 29) |
			(vdp->mantissa3 << 13) |
			((vdp->mantissa4 >> 3) & MASK(13));

		unsigned mant_hi = (vdp->mantissa1 << 13)
				 | (vdp->mantissa2 >> 3);

		if((vdp->mantissa4 & 7) > 4)
		{
			/* round up */
			mant_lo++;
			if(mant_lo == 0)
			{
				mant_hi++;
				if(mant_hi > 0xffffff)
				{
					mant_hi = 0;
					exp++;
				}
			}
		}

		idp->mant_lo = SWAP4(mant_lo);
		idp->mant_6 = mant_hi >> 16;
		idp->mant_5 = (mant_hi & 0xff00) >> 8;
		idp->mant_4 = mant_hi;
		idp->exp_hi = exp >> 4;
		idp->exp_lo = exp;
	}
		
	shipit:
		idp->sign = vdp->sign;

}

	/* vax */
#else
#error "ix_double implementation"
#endif

int
ncmpix_get_double_schar(const void *xp, schar *ip)
{
	double xx=0.0;
	get_ix_double(xp, &xx);
	*ip = (schar) xx;
	if(xx > SCHAR_MAX || xx < SCHAR_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_double_uchar(const void *xp, uchar *ip)
{
	double xx=0.0;
	get_ix_double(xp, &xx);
	*ip = (uchar) xx;
	if(xx > UCHAR_MAX || xx < 0)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_double_short(const void *xp, short *ip)
{
	double xx=0.0;
	get_ix_double(xp, &xx);
	*ip = (short) xx;
	if(xx > SHORT_MAX || xx < SHORT_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_double_int(const void *xp, int *ip)
{
	double xx=0.0;
	get_ix_double(xp, &xx);
	*ip = (int) xx;
	if(xx > INT_MAX || xx < INT_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_double_long(const void *xp, long *ip)
{
	double xx;
	get_ix_double(xp, &xx);
	*ip = (long) xx;
	if(xx > LONG_MAX || xx < LONG_MIN)
		return NC_ERANGE;
	return NC_NOERR;
}

int
ncmpix_get_double_float(const void *xp, float *ip)
{
	double xx=0.0;
	get_ix_double(xp, &xx);
	if(xx > FLT_MAX || xx < (-FLT_MAX))
	{
		*ip = FLT_MAX;
		return NC_ERANGE;
	}
	if(xx < (-FLT_MAX))
	{
		*ip = (-FLT_MAX);
		return NC_ERANGE;
	}
	*ip = (float) xx;
	return NC_NOERR;
}

int
ncmpix_get_double_double(const void *xp, double *ip)
{
	/* TODO */
	get_ix_double(xp, ip);
	return NC_NOERR;
}


int
ncmpix_put_double_schar(void *xp, const schar *ip)
{
	double xx = (double) *ip;
	put_ix_double(xp, &xx);
	return NC_NOERR;
}

int
ncmpix_put_double_uchar(void *xp, const uchar *ip)
{
	double xx = (double) *ip;
	put_ix_double(xp, &xx);
	return NC_NOERR;
}

int
ncmpix_put_double_short(void *xp, const short *ip)
{
	double xx = (double) *ip;
	put_ix_double(xp, &xx);
#if 0	/* TODO: figure this out */
	if((double)(*ip) > X_DOUBLE_MAX || (double)(*ip) < X_DOUBLE_MIN)
		return NC_ERANGE;
#endif
	return NC_NOERR;
}

int
ncmpix_put_double_int(void *xp, const int *ip)
{
	double xx = (double) *ip;
	put_ix_double(xp, &xx);
#if 0	/* TODO: figure this out */
	if((double)(*ip) > X_DOUBLE_MAX || (double)(*ip) < X_DOUBLE_MIN)
		return NC_ERANGE;
#endif
	return NC_NOERR;
}

int
ncmpix_put_double_long(void *xp, const long *ip)
{
	double xx = (double) *ip;
	put_ix_double(xp, &xx);
#if 1	/* TODO: figure this out */
	if((double)(*ip) > X_DOUBLE_MAX || (double)(*ip) < X_DOUBLE_MIN)
		return NC_ERANGE;
#endif
	return NC_NOERR;
}

int
ncmpix_put_double_float(void *xp, const float *ip)
{
	double xx = (double) *ip;
	put_ix_double(xp, &xx);
#if 1	/* TODO: figure this out */
	if((double)(*ip) > X_DOUBLE_MAX || (double)(*ip) < X_DOUBLE_MIN)
		return NC_ERANGE;
#endif
	return NC_NOERR;
}

int
ncmpix_put_double_double(void *xp, const double *ip)
{
	put_ix_double(xp, ip);
#ifdef NO_IEEE_FLOAT
	if(*ip > X_DOUBLE_MAX || *ip < X_DOUBLE_MIN)
		return NC_ERANGE;
#endif
	return NC_NOERR;
}


/* x_size_t */

#if SIZEOF_SIZE_T < X_SIZEOF_SIZE_T
#error "x_size_t implementation"
/* netcdf requires MPI_Offset which can hold a values from 0 to 2^31 -1 */
#endif

int
ncmpix_put_size_t1(void **xpp, const MPI_Offset *ulp)
{
	/* similar to put_ix_int() */
	uchar *cp = (uchar *) *xpp;
	assert(*ulp <= X_SIZE_MAX);

	*cp++ = (uchar)((*ulp) >> 24);
	*cp++ = (uchar)(((*ulp) & 0x00ff0000) >> 16);
	*cp++ = (uchar)(((*ulp) & 0x0000ff00) >>  8);
	*cp   = (uchar)((*ulp) & 0x000000ff);

	*xpp = (void *)((char *)(*xpp) + X_SIZEOF_SIZE_T);
	return NC_NOERR;
}

/*
int
ncmpix_put_size_t(void **xpp, const MPI_Offset *ulp)
{
//	uchar *cp = (uchar *) *xpp;
	unsigned long *cp = (unsigned long *) *xpp;
	assert(*ulp <= X_SIZE_MAX);

	*cp++ = (uchar)((*ulp) >> 24);
	*cp++ = (uchar)(((*ulp) & 0x00ff0000) >> 16);
	*cp++ = (uchar)(((*ulp) & 0x0000ff00) >>  8);
	*cp   = (uchar)((*ulp) & 0x000000ff);

	*xpp = (void *)((char *)(*xpp) + X_SIZEOF_LONG);
	return NC_NOERR;
}
*/
int
ncmpix_put_size_t(void **xpp, const MPI_Offset *lp, MPI_Offset sizeof_t)
{
        /* similar to put_ix_int() */
        uchar *cp = (uchar *) *xpp;
//	assert(*lp <= X_SIZE_MAX);
                /* No negative offsets stored in netcdf */
//        if (*lp < 0) {
                /* assume this is an overflow of a 32-bit int */
//               return ERANGE;
//        }

        assert(sizeof_t == 4 || sizeof_t == 8);
        if (sizeof_t == 4 ) {
                *cp++ = (uchar)(((*lp) & 0xff000000) >> 24);
                *cp++ = (uchar)(((*lp) & 0x00ff0000) >> 16);
                *cp++ = (uchar)(((*lp) & 0x0000ff00) >>  8);
                *cp   = (uchar)( (*lp) & 0x000000ff);
	        *xpp = (void *)((char *)(*xpp) + sizeof_t);
//	        *xpp = (void *)((char *)(*xpp) + X_SIZEOF_SIZE_T);
        } else {
                *cp++ = (uchar)(((*lp) & 0xff00000000000000ULL) >> 56);
                *cp++ = (uchar)(((*lp) & 0x00ff000000000000ULL) >> 48);
                *cp++ = (uchar)(((*lp) & 0x0000ff0000000000ULL) >> 40);
                *cp++ = (uchar)(((*lp) & 0x000000ff00000000ULL) >> 32);
                *cp++ = (uchar)(((*lp) & 0x00000000ff000000ULL) >> 24);
                *cp++ = (uchar)(((*lp) & 0x0000000000ff0000ULL) >> 16);
                *cp++ = (uchar)(((*lp) & 0x000000000000ff00ULL) >>  8);
                *cp   = (uchar)( (*lp) & 0x00000000000000ffULL);
	        *xpp = (void *)((char *)(*xpp) + sizeof_t);
        }
        return NC_NOERR;
}




int
ncmpix_get_size_t(const void **xpp,  MPI_Offset *lp, MPI_Offset sizeof_off_t)
{
	/* similar to get_ix_int */
	/* similar to get_ix_int() */
	const uchar *cp = (const uchar *) *xpp;
	if (*lp < 0) {
		/* assume this is an overflow of a 32-bit int */
		return ERANGE;
	}
	assert(sizeof_off_t == 4 || sizeof_off_t == 8);
       if (sizeof_off_t == 4) {
               *lp = *cp++ << 24;
               *lp |= (*cp++ << 16);
               *lp |= (*cp++ <<  8);
               *lp |= *cp; 
       } else {
               *lp =  ((off_t)(*cp++) << 56);
               *lp |= ((off_t)(*cp++) << 48);
               *lp |= ((off_t)(*cp++) << 40);
               *lp |= ((off_t)(*cp++) << 32);
               *lp |= ((off_t)(*cp++) << 24);
               *lp |= ((off_t)(*cp++) << 16);
               *lp |= ((off_t)(*cp++) <<  8);
               *lp |=  (off_t)*cp;
       }
	*xpp = (const void *)((const char *)(*xpp) + sizeof_off_t);
	return NC_NOERR;
}

/* x_off_t */

/* A previous version of this function would try to let systems with a 4 byte
 * off_t read and write the 8-byte offsets of a CDF-2 file.  For simplicity, we
 * now ensure that platforms with a 4-byte off_t will never open a CDF-2 file
 * no matter what the size. 
 */
int
ncmpix_put_off_t(void **xpp, const MPI_Offset *lp, MPI_Offset sizeof_off_t)
{
	/* similar to put_ix_int() */
	uchar *cp = (uchar *) *xpp;
		/* No negative offsets stored in netcdf */
	if (*lp < 0) {
		/* assume this is an overflow of a 32-bit int */
		return ERANGE;
	}
	assert(sizeof_off_t == 4 || sizeof_off_t == 8);
	if (sizeof_off_t == 4 ) {
		*cp++ = (uchar)(((*lp) & 0xff000000) >> 24);
		*cp++ = (uchar)(((*lp) & 0x00ff0000) >> 16);
		*cp++ = (uchar)(((*lp) & 0x0000ff00) >>  8);
		*cp   = (uchar)( (*lp) & 0x000000ff);
	} else {
		*cp++ = (uchar)(((*lp) & 0xff00000000000000ULL) >> 56);
		*cp++ = (uchar)(((*lp) & 0x00ff000000000000ULL) >> 48);
		*cp++ = (uchar)(((*lp) & 0x0000ff0000000000ULL) >> 40);
		*cp++ = (uchar)(((*lp) & 0x000000ff00000000ULL) >> 32);
		*cp++ = (uchar)(((*lp) & 0x00000000ff000000ULL) >> 24);
		*cp++ = (uchar)(((*lp) & 0x0000000000ff0000ULL) >> 16);
		*cp++ = (uchar)(((*lp) & 0x000000000000ff00ULL) >>  8);
		*cp   = (uchar)( (*lp) & 0x00000000000000ffULL);
	}
	*xpp = (void *)((char *)(*xpp) + sizeof_off_t);
	return NC_NOERR;
}

/* see comments for ncmpix_put_off_t */
int
ncmpix_get_off_t(const void **xpp, MPI_Offset *lp, MPI_Offset sizeof_off_t)
{
	/* similar to get_ix_int() */
	const uchar *cp = (const uchar *) *xpp;
	assert(sizeof_off_t == 4 || sizeof_off_t == 8);
	
       if (sizeof_off_t == 4) {
               *lp = *cp++ << 24;
               *lp |= (*cp++ << 16);
               *lp |= (*cp++ <<  8);
               *lp |= *cp; 
       } else {
               *lp =  ((off_t)(*cp++) << 56);
               *lp |= ((off_t)(*cp++) << 48);
               *lp |= ((off_t)(*cp++) << 40);
               *lp |= ((off_t)(*cp++) << 32);
               *lp |= ((off_t)(*cp++) << 24);
               *lp |= ((off_t)(*cp++) << 16);
               *lp |= ((off_t)(*cp++) <<  8);
               *lp |=  (off_t)*cp;
       }
       *xpp = (const void *)((const char *)(*xpp) + sizeof_off_t);


	return NC_NOERR;
}


/*
 * Aggregate numeric conversion functions.
 */



/* schar */

int
ncmpix_getn_schar_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
		(void) memcpy(tp, *xpp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);
	return NC_NOERR;

}
int
ncmpix_getn_schar_uchar(const void **xpp, MPI_Offset nelems, uchar *tp)
{
		(void) memcpy(tp, *xpp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);
	return NC_NOERR;

}
int
ncmpix_getn_schar_short(const void **xpp, MPI_Offset nelems, short *tp)
{
	schar *xp = (schar *)(*xpp);

	while(nelems-- != 0)
	{
		*tp++ = *xp++;
	}

	*xpp = (const void *)xp;
	return NC_NOERR;
}

int
ncmpix_getn_schar_int(const void **xpp, MPI_Offset nelems, int *tp)
{
	schar *xp = (schar *)(*xpp);

	while(nelems-- != 0)
	{
		*tp++ = *xp++;
	}

	*xpp = (const void *)xp;
	return NC_NOERR;
}

int
ncmpix_getn_schar_long(const void **xpp, MPI_Offset nelems, long *tp)
{
	schar *xp = (schar *)(*xpp);

	while(nelems-- != 0)
	{
		*tp++ = *xp++;
	}

	*xpp = (const void *)xp;
	return NC_NOERR;
}

int
ncmpix_getn_schar_float(const void **xpp, MPI_Offset nelems, float *tp)
{
	schar *xp = (schar *)(*xpp);

	while(nelems-- != 0)
	{
		*tp++ = *xp++;
	}

	*xpp = (const void *)xp;
	return NC_NOERR;
}

int
ncmpix_getn_schar_double(const void **xpp, MPI_Offset nelems, double *tp)
{
	schar *xp = (schar *)(*xpp);

	while(nelems-- != 0)
	{
		*tp++ = *xp++;
	}

	*xpp = (const void *)xp;
	return NC_NOERR;
}


int
ncmpix_pad_getn_schar_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
		MPI_Offset rndup = nelems % X_ALIGN;

	if(rndup)
		rndup = X_ALIGN - rndup;

	(void) memcpy(tp, *xpp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems + rndup);

	return NC_NOERR;

}
int
ncmpix_pad_getn_schar_uchar(const void **xpp, MPI_Offset nelems, uchar *tp)
{
		MPI_Offset rndup = nelems % X_ALIGN;

	if(rndup)
		rndup = X_ALIGN - rndup;

	(void) memcpy(tp, *xpp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems + rndup);

	return NC_NOERR;

}
int
ncmpix_pad_getn_schar_short(const void **xpp, MPI_Offset nelems, short *tp)
{
	MPI_Offset rndup = nelems % X_ALIGN;
	schar *xp = (schar *) *xpp;

	if(rndup)
		rndup = X_ALIGN - rndup;

	while(nelems-- != 0)
	{
		*tp++ = *xp++;
	}

	*xpp = (void *)(xp + rndup);
	return NC_NOERR;
}

int
ncmpix_pad_getn_schar_int(const void **xpp, MPI_Offset nelems, int *tp)
{
	MPI_Offset rndup = nelems % X_ALIGN;
	schar *xp = (schar *) *xpp;

	if(rndup)
		rndup = X_ALIGN - rndup;

	while(nelems-- != 0)
	{
		*tp++ = *xp++;
	}

	*xpp = (void *)(xp + rndup);
	return NC_NOERR;
}

int
ncmpix_pad_getn_schar_long(const void **xpp, MPI_Offset nelems, long *tp)
{
	MPI_Offset rndup = nelems % X_ALIGN;
	schar *xp = (schar *) *xpp;

	if(rndup)
		rndup = X_ALIGN - rndup;

	while(nelems-- != 0)
	{
		*tp++ = *xp++;
	}

	*xpp = (void *)(xp + rndup);
	return NC_NOERR;
}

int
ncmpix_pad_getn_schar_float(const void **xpp, MPI_Offset nelems, float *tp)
{
	MPI_Offset rndup = nelems % X_ALIGN;
	schar *xp = (schar *) *xpp;

	if(rndup)
		rndup = X_ALIGN - rndup;

	while(nelems-- != 0)
	{
		*tp++ = *xp++;
	}

	*xpp = (void *)(xp + rndup);
	return NC_NOERR;
}

int
ncmpix_pad_getn_schar_double(const void **xpp, MPI_Offset nelems, double *tp)
{
	MPI_Offset rndup = nelems % X_ALIGN;
	schar *xp = (schar *) *xpp;

	if(rndup)
		rndup = X_ALIGN - rndup;

	while(nelems-- != 0)
	{
		*tp++ = *xp++;
	}

	*xpp = (void *)(xp + rndup);
	return NC_NOERR;
}


int
ncmpix_putn_schar_schar(void **xpp, MPI_Offset nelems, const schar *tp)
{
		(void) memcpy(*xpp, tp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);

	return NC_NOERR;

}
int
ncmpix_putn_schar_uchar(void **xpp, MPI_Offset nelems, const uchar *tp)
{
		(void) memcpy(*xpp, tp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);

	return NC_NOERR;

}
int
ncmpix_putn_schar_short(void **xpp, MPI_Offset nelems, const short *tp)
{
	int status = NC_NOERR;
	schar *xp = (schar *) *xpp;

	while(nelems-- != 0)
	{
		if(*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN)
			status = NC_ERANGE;
		*xp++ = (schar) *tp++;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_schar_int(void **xpp, MPI_Offset nelems, const int *tp)
{
	int status = NC_NOERR;
	schar *xp = (schar *) *xpp;

	while(nelems-- != 0)
	{
		if(*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN)
			status = NC_ERANGE;
		*xp++ = (schar) *tp++;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_schar_long(void **xpp, MPI_Offset nelems, const long *tp)
{
	int status = NC_NOERR;
	schar *xp = (schar *) *xpp;

	while(nelems-- != 0)
	{
		if(*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN)
			status = NC_ERANGE;
		*xp++ = (schar) *tp++;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_schar_float(void **xpp, MPI_Offset nelems, const float *tp)
{
	int status = NC_NOERR;
	schar *xp = (schar *) *xpp;

	while(nelems-- != 0)
	{
		if(*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN)
			status = NC_ERANGE;
		*xp++ = (schar) *tp++;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_schar_double(void **xpp, MPI_Offset nelems, const double *tp)
{
	int status = NC_NOERR;
	schar *xp = (schar *) *xpp;

	while(nelems-- != 0)
	{
		if(*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN)
			status = NC_ERANGE;
		*xp++ = (schar) *tp++;
	}

	*xpp = (void *)xp;
	return status;
}


int
ncmpix_pad_putn_schar_schar(void **xpp, MPI_Offset nelems, const schar *tp)
{
		MPI_Offset rndup = nelems % X_ALIGN;

	if(rndup)
		rndup = X_ALIGN - rndup;

	(void) memcpy(*xpp, tp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);

	if(rndup)
	{
		(void) memcpy(*xpp, nada, rndup);
		*xpp = (void *)((char *)(*xpp) + rndup);
	}
	
	return NC_NOERR;

}
int
ncmpix_pad_putn_schar_uchar(void **xpp, MPI_Offset nelems, const uchar *tp)
{
		MPI_Offset rndup = nelems % X_ALIGN;

	if(rndup)
		rndup = X_ALIGN - rndup;

	(void) memcpy(*xpp, tp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);

	if(rndup)
	{
		(void) memcpy(*xpp, nada, rndup);
		*xpp = (void *)((char *)(*xpp) + rndup);
	}
	
	return NC_NOERR;

}
int
ncmpix_pad_putn_schar_short(void **xpp, MPI_Offset nelems, const short *tp)
{
	int status = NC_NOERR;
	MPI_Offset rndup = nelems % X_ALIGN;
	schar *xp = (schar *) *xpp;

	if(rndup)
		rndup = X_ALIGN - rndup;

	while(nelems-- != 0)
	{
		/* N.B. schar as signed */
		if(*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN)
			status = NC_ERANGE;
		*xp++ = (schar) *tp++;
	}


	if(rndup)
	{
		(void) memcpy(xp, nada, rndup);
		xp += rndup;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_putn_schar_int(void **xpp, MPI_Offset nelems, const int *tp)
{
	int status = NC_NOERR;
	MPI_Offset rndup = nelems % X_ALIGN;
	schar *xp = (schar *) *xpp;

	if(rndup)
		rndup = X_ALIGN - rndup;

	while(nelems-- != 0)
	{
		/* N.B. schar as signed */
		if(*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN)
			status = NC_ERANGE;
		*xp++ = (schar) *tp++;
	}


	if(rndup)
	{
		(void) memcpy(xp, nada, rndup);
		xp += rndup;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_putn_schar_long(void **xpp, MPI_Offset nelems, const long *tp)
{
	int status = NC_NOERR;
	MPI_Offset rndup = nelems % X_ALIGN;
	schar *xp = (schar *) *xpp;

	if(rndup)
		rndup = X_ALIGN - rndup;

	while(nelems-- != 0)
	{
		/* N.B. schar as signed */
		if(*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN)
			status = NC_ERANGE;
		*xp++ = (schar) *tp++;
	}


	if(rndup)
	{
		(void) memcpy(xp, nada, rndup);
		xp += rndup;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_putn_schar_float(void **xpp, MPI_Offset nelems, const float *tp)
{
	int status = NC_NOERR;
	MPI_Offset rndup = nelems % X_ALIGN;
	schar *xp = (schar *) *xpp;

	if(rndup)
		rndup = X_ALIGN - rndup;

	while(nelems-- != 0)
	{
		/* N.B. schar as signed */
		if(*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN)
			status = NC_ERANGE;
		*xp++ = (schar) *tp++;
	}


	if(rndup)
	{
		(void) memcpy(xp, nada, rndup);
		xp += rndup;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_putn_schar_double(void **xpp, MPI_Offset nelems, const double *tp)
{
	int status = NC_NOERR;
	MPI_Offset rndup = nelems % X_ALIGN;
	schar *xp = (schar *) *xpp;

	if(rndup)
		rndup = X_ALIGN - rndup;

	while(nelems-- != 0)
	{
		/* N.B. schar as signed */
		if(*tp > X_SCHAR_MAX || *tp < X_SCHAR_MIN)
			status = NC_ERANGE;
		*xp++ = (schar) *tp++;
	}


	if(rndup)
	{
		(void) memcpy(xp, nada, rndup);
		xp += rndup;
	}

	*xpp = (void *)xp;
	return status;
}



/* short */

int
ncmpix_getn_short_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_schar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_short_uchar(const void **xpp, MPI_Offset nelems, uchar *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_uchar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

#if X_SIZEOF_SHORT == SIZEOF_SHORT
/* optimized version */
int
ncmpix_getn_short_short(const void **xpp, MPI_Offset nelems, short *tp)
{
# ifdef WORDS_BIGENDIAN
	(void) memcpy(tp, *xpp, nelems * sizeof(short));
# else
	swapn2b(tp, *xpp, nelems);
# endif
	*xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_SHORT);
	return NC_NOERR;
}
#else
int
ncmpix_getn_short_short(const void **xpp, MPI_Offset nelems, short *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_short(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

#endif
int
ncmpix_getn_short_int(const void **xpp, MPI_Offset nelems, int *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_int(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_short_long(const void **xpp, MPI_Offset nelems, long *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_long(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_short_float(const void **xpp, MPI_Offset nelems, float *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_float(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_short_double(const void **xpp, MPI_Offset nelems, double *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_double(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}


int
ncmpix_pad_getn_short_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
	const MPI_Offset rndup = nelems % 2;

	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_schar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
		xp += X_SIZEOF_SHORT;
		
	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_getn_short_uchar(const void **xpp, MPI_Offset nelems, uchar *tp)
{
	const MPI_Offset rndup = nelems % 2;

	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_uchar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
		xp += X_SIZEOF_SHORT;
		
	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_getn_short_short(const void **xpp, MPI_Offset nelems, short *tp)
{
	const MPI_Offset rndup = nelems % 2;

	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_short(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
		xp += X_SIZEOF_SHORT;
		
	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_getn_short_int(const void **xpp, MPI_Offset nelems, int *tp)
{
	const MPI_Offset rndup = nelems % 2;

	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_int(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
		xp += X_SIZEOF_SHORT;
		
	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_getn_short_long(const void **xpp, MPI_Offset nelems, long *tp)
{
	const MPI_Offset rndup = nelems % 2;

	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_long(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
		xp += X_SIZEOF_SHORT;
		
	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_getn_short_float(const void **xpp, MPI_Offset nelems, float *tp)
{
	const MPI_Offset rndup = nelems % 2;

	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_float(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
		xp += X_SIZEOF_SHORT;
		
	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_getn_short_double(const void **xpp, MPI_Offset nelems, double *tp)
{
	const MPI_Offset rndup = nelems % 2;

	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		const int lstatus = ncmpix_get_short_double(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
		xp += X_SIZEOF_SHORT;
		
	*xpp = (void *)xp;
	return status;
}


int
ncmpix_putn_short_schar(void **xpp, MPI_Offset nelems, const schar *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_schar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_short_uchar(void **xpp, MPI_Offset nelems, const uchar *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_uchar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

#if X_SIZEOF_SHORT == SIZEOF_SHORT
/* optimized version */
int
ncmpix_putn_short_short(void **xpp, MPI_Offset nelems, const short *tp)
{
# ifdef WORDS_BIGENDIAN
	(void) memcpy(*xpp, tp, nelems * X_SIZEOF_SHORT);
# else
	swapn2b(*xpp, tp, nelems);
# endif
	*xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_SHORT);
	return NC_NOERR;
}
#else
int
ncmpix_putn_short_short(void **xpp, MPI_Offset nelems, const short *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_short(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

#endif
int
ncmpix_putn_short_int(void **xpp, MPI_Offset nelems, const int *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_int(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_short_long(void **xpp, MPI_Offset nelems, const long *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_long(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_short_float(void **xpp, MPI_Offset nelems, const float *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_float(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_short_double(void **xpp, MPI_Offset nelems, const double *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_double(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}


int
ncmpix_pad_putn_short_schar(void **xpp, MPI_Offset nelems, const schar *tp)
{
	const MPI_Offset rndup = nelems % 2;

	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_schar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
	{
		(void) memcpy(xp, nada, X_SIZEOF_SHORT);
		xp += X_SIZEOF_SHORT;	
	}
		
	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_putn_short_uchar(void **xpp, MPI_Offset nelems, const uchar *tp)
{
	const MPI_Offset rndup = nelems % 2;

	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_uchar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
	{
		(void) memcpy(xp, nada, X_SIZEOF_SHORT);
		xp += X_SIZEOF_SHORT;	
	}
		
	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_putn_short_short(void **xpp, MPI_Offset nelems, const short *tp)
{
	const MPI_Offset rndup = nelems % 2;

	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_short(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
	{
		(void) memcpy(xp, nada, X_SIZEOF_SHORT);
		xp += X_SIZEOF_SHORT;	
	}
		
	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_putn_short_int(void **xpp, MPI_Offset nelems, const int *tp)
{
	const MPI_Offset rndup = nelems % 2;

	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_int(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
	{
		(void) memcpy(xp, nada, X_SIZEOF_SHORT);
		xp += X_SIZEOF_SHORT;	
	}
		
	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_putn_short_long(void **xpp, MPI_Offset nelems, const long *tp)
{
	const MPI_Offset rndup = nelems % 2;

	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_long(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
	{
		(void) memcpy(xp, nada, X_SIZEOF_SHORT);
		xp += X_SIZEOF_SHORT;	
	}
		
	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_putn_short_float(void **xpp, MPI_Offset nelems, const float *tp)
{
	const MPI_Offset rndup = nelems % 2;

	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_float(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
	{
		(void) memcpy(xp, nada, X_SIZEOF_SHORT);
		xp += X_SIZEOF_SHORT;	
	}
		
	*xpp = (void *)xp;
	return status;
}

int
ncmpix_pad_putn_short_double(void **xpp, MPI_Offset nelems, const double *tp)
{
	const MPI_Offset rndup = nelems % 2;

	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_SHORT, tp++)
	{
		int lstatus = ncmpix_put_short_double(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	if(rndup != 0)
	{
		(void) memcpy(xp, nada, X_SIZEOF_SHORT);
		xp += X_SIZEOF_SHORT;	
	}
		
	*xpp = (void *)xp;
	return status;
}



/* int */

int
ncmpix_getn_int_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		const int lstatus = ncmpix_get_int_schar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_int_uchar(const void **xpp, MPI_Offset nelems, uchar *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		const int lstatus = ncmpix_get_int_uchar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_int_short(const void **xpp, MPI_Offset nelems, short *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		const int lstatus = ncmpix_get_int_short(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

#if X_SIZEOF_INT == SIZEOF_INT
/* optimized version */
int
ncmpix_getn_int_int(const void **xpp, MPI_Offset nelems, int *tp)
{
# ifdef WORDS_BIGENDIAN
	(void) memcpy(tp, *xpp, nelems * sizeof(int));
# else
	swapn4b(tp, *xpp, nelems);
# endif
	*xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_INT);
	return NC_NOERR;
}
#else
int
ncmpix_getn_int_int(const void **xpp, MPI_Offset nelems, int *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		const int lstatus = ncmpix_get_int_int(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

#endif
#if X_SIZEOF_INT == SIZEOF_LONG
/* optimized version */
int
ncmpix_getn_int_long(const void **xpp, MPI_Offset nelems, long *tp)
{
# ifdef WORDS_BIGENDIAN
	(void) memcpy(tp, *xpp, nelems * sizeof(long));
# else
	swapn4b(tp, *xpp, nelems);
# endif
	*xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_INT);
	return NC_NOERR;
}
#else
int
ncmpix_getn_int_long(const void **xpp, MPI_Offset nelems, long *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		const int lstatus = ncmpix_get_int_long(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

#endif

int
ncmpix_getn_long_long(const void **xpp, MPI_Offset nelems, MPI_Offset *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += 8, tp++)
	{
		const int lstatus = ncmpix_get_long_long(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}


int
ncmpix_getn_int_float(const void **xpp, MPI_Offset nelems, float *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		const int lstatus = ncmpix_get_int_float(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_int_double(const void **xpp, MPI_Offset nelems, double *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		const int lstatus = ncmpix_get_int_double(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}


int
ncmpix_putn_int_schar(void **xpp, MPI_Offset nelems, const schar *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		int lstatus = ncmpix_put_int_schar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_int_uchar(void **xpp, MPI_Offset nelems, const uchar *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		int lstatus = ncmpix_put_int_uchar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_int_short(void **xpp, MPI_Offset nelems, const short *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		int lstatus = ncmpix_put_int_short(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

#if X_SIZEOF_INT == SIZEOF_INT
/* optimized version */
int
ncmpix_putn_int_int(void **xpp, MPI_Offset nelems, const int *tp)
{
# ifdef WORDS_BIGENDIAN
	(void) memcpy(*xpp, tp, nelems * X_SIZEOF_INT);
# else
	swapn4b(*xpp, tp, nelems);
# endif
	*xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_INT);
	return NC_NOERR;
}
#else
int
ncmpix_putn_int_int(void **xpp, MPI_Offset nelems, const int *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		int lstatus = ncmpix_put_int_int(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

#endif
#if X_SIZEOF_INT == SIZEOF_LONG
/* optimized version */
int
ncmpix_putn_int_long(void **xpp, MPI_Offset nelems, const long *tp)
{
# ifdef WORDS_BIGENDIAN
	(void) memcpy(*xpp, tp, nelems * X_SIZEOF_INT);
# else
	swapn4b(*xpp, tp, nelems);
# endif
	*xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_INT);
	return NC_NOERR;
}
#else
int
ncmpix_putn_int_long(void **xpp, MPI_Offset nelems, const long *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		int lstatus = ncmpix_put_int_long(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

#endif
int
ncmpix_putn_int_float(void **xpp, MPI_Offset nelems, const float *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		int lstatus = ncmpix_put_int_float(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_int_double(void **xpp, MPI_Offset nelems, const double *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_INT, tp++)
	{
		int lstatus = ncmpix_put_int_double(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}



/* float */

int
ncmpix_getn_float_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		const int lstatus = ncmpix_get_float_schar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_float_uchar(const void **xpp, MPI_Offset nelems, uchar *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		const int lstatus = ncmpix_get_float_uchar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_float_short(const void **xpp, MPI_Offset nelems, short *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		const int lstatus = ncmpix_get_float_short(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_float_int(const void **xpp, MPI_Offset nelems, int *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		const int lstatus = ncmpix_get_float_int(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_float_long(const void **xpp, MPI_Offset nelems, long *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		const int lstatus = ncmpix_get_float_long(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

#if X_SIZEOF_FLOAT == SIZEOF_FLOAT && !defined(NO_IEEE_FLOAT)
/* optimized version */
int
ncmpix_getn_float_float(const void **xpp, MPI_Offset nelems, float *tp)
{
# ifdef WORDS_BIGENDIAN
	(void) memcpy(tp, *xpp, nelems * sizeof(float));
# else
	swapn4b(tp, *xpp, nelems);
# endif
	*xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_FLOAT);
	return NC_NOERR;
}
#elif vax
int
ncmpix_getn_float_float(const void **xpp, MPI_Offset nfloats, float *ip)
{
	float *const end = ip + nfloats;

	while(ip < end)
	{
		struct vax_single *const vsp = (struct vax_single *) ip;
		const struct ieee_single *const isp =
			 (const struct ieee_single *) (*xpp);
		unsigned exp = isp->exp_hi << 1 | isp->exp_lo;

		switch(exp) {
		case 0 :
			/* ieee subnormal */
			if(isp->mant_hi == min.ieee.mant_hi
				&& isp->mant_lo_hi == min.ieee.mant_lo_hi
				&& isp->mant_lo_lo == min.ieee.mant_lo_lo)
			{
				*vsp = min.s;
			}
			else
			{
				unsigned mantissa = (isp->mant_hi << 16)
					 | isp->mant_lo_hi << 8
					 | isp->mant_lo_lo;
				unsigned tmp = mantissa >> 20;
				if(tmp >= 4) {
					vsp->exp = 2;
				} else if (tmp >= 2) {
					vsp->exp = 1;
				} else {
					*vsp = min.s;
					break;
				} /* else */
				tmp = mantissa - (1 << (20 + vsp->exp ));
				tmp <<= 3 - vsp->exp;
				vsp->mantissa2 = tmp;
				vsp->mantissa1 = (tmp >> 16);
			}
			break;
		case 0xfe :
		case 0xff :
			*vsp = max.s;
			break;
		default :
			vsp->exp = exp - IEEE_SNG_BIAS + VAX_SNG_BIAS;
			vsp->mantissa2 = isp->mant_lo_hi << 8 | isp->mant_lo_lo;
			vsp->mantissa1 = isp->mant_hi;
		}

		vsp->sign = isp->sign;


		ip++;
		*xpp = (char *)(*xpp) + X_SIZEOF_FLOAT;
	}
	return NC_NOERR;
}
#else
int
ncmpix_getn_float_float(const void **xpp, MPI_Offset nelems, float *tp)
{
	const char *xp = *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		const int lstatus = ncmpix_get_float_float(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

#endif
int
ncmpix_getn_float_double(const void **xpp, MPI_Offset nelems, double *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		const int lstatus = ncmpix_get_float_double(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}


int
ncmpix_putn_float_schar(void **xpp, MPI_Offset nelems, const schar *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		int lstatus = ncmpix_put_float_schar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_float_uchar(void **xpp, MPI_Offset nelems, const uchar *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		int lstatus = ncmpix_put_float_uchar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_float_short(void **xpp, MPI_Offset nelems, const short *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		int lstatus = ncmpix_put_float_short(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_float_int(void **xpp, MPI_Offset nelems, const int *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		int lstatus = ncmpix_put_float_int(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_float_long(void **xpp, MPI_Offset nelems, const long *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		int lstatus = ncmpix_put_float_long(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

#if X_SIZEOF_FLOAT == SIZEOF_FLOAT && !defined(NO_IEEE_FLOAT)
/* optimized version */
int
ncmpix_putn_float_float(void **xpp, MPI_Offset nelems, const float *tp)
{
# ifdef WORDS_BIGENDIAN
	(void) memcpy(*xpp, tp, nelems * X_SIZEOF_FLOAT);
# else
	swapn4b(*xpp, tp, nelems);
# endif
	*xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_FLOAT);
	return NC_NOERR;
}
#elif vax
int
ncmpix_putn_float_float(void **xpp, MPI_Offset nfloats, const float *ip)
{
	const float *const end = ip + nfloats;

	while(ip < end)
	{
		const struct vax_single *const vsp =
			 (const struct vax_single *)ip;
		struct ieee_single *const isp = (struct ieee_single *) (*xpp);

		switch(vsp->exp){
		case 0 :
			/* all vax float with zero exponent map to zero */
			*isp = min.ieee;
			break;
		case 2 :
		case 1 :
		{
			/* These will map to subnormals */
			unsigned mantissa = (vsp->mantissa1 << 16)
					 | vsp->mantissa2;
			mantissa >>= 3 - vsp->exp;
			mantissa += (1 << (20 + vsp->exp));
			isp->mant_lo_lo = mantissa;
			isp->mant_lo_hi = mantissa >> 8;
			isp->mant_hi = mantissa >> 16;
			isp->exp_lo = 0;
			isp->exp_hi = 0;
		}
			break;
		case 0xff : /* max.s.exp */
			if( vsp->mantissa2 == max.s.mantissa2
				&& vsp->mantissa1 == max.s.mantissa1)
			{
				/* map largest vax float to ieee infinity */
				*isp = max.ieee;
				break;
			} /* else, fall thru */
		default :
		{
			unsigned exp = vsp->exp - VAX_SNG_BIAS + IEEE_SNG_BIAS;
			isp->exp_hi = exp >> 1;
			isp->exp_lo = exp;
			isp->mant_lo_lo = vsp->mantissa2;
			isp->mant_lo_hi = vsp->mantissa2 >> 8;
			isp->mant_hi = vsp->mantissa1;
		}
		}

		isp->sign = vsp->sign;

	
		ip++;
		*xpp = (char *)(*xpp) + X_SIZEOF_FLOAT;
	}
	return NC_NOERR;
}
#else
int
ncmpix_putn_float_float(void **xpp, MPI_Offset nelems, const float *tp)
{
	char *xp = *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		int lstatus = ncmpix_put_float_float(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

#endif
int
ncmpix_putn_float_double(void **xpp, MPI_Offset nelems, const double *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++)
	{
		int lstatus = ncmpix_put_float_double(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}



/* double */

int
ncmpix_getn_double_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		const int lstatus = ncmpix_get_double_schar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_double_uchar(const void **xpp, MPI_Offset nelems, uchar *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		const int lstatus = ncmpix_get_double_uchar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_double_short(const void **xpp, MPI_Offset nelems, short *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		const int lstatus = ncmpix_get_double_short(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_double_int(const void **xpp, MPI_Offset nelems, int *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		const int lstatus = ncmpix_get_double_int(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_double_long(const void **xpp, MPI_Offset nelems, long *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		const int lstatus = ncmpix_get_double_long(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

int
ncmpix_getn_double_float(const void **xpp, MPI_Offset nelems, float *tp)
{
	const char *xp = (const char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		const int lstatus = ncmpix_get_double_float(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

#if X_SIZEOF_DOUBLE == SIZEOF_DOUBLE && !defined(NO_IEEE_FLOAT)
/* optimized version */
int
ncmpix_getn_double_double(const void **xpp, MPI_Offset nelems, double *tp)
{
# ifdef WORDS_BIGENDIAN
	(void) memcpy(tp, *xpp, nelems * sizeof(double));
# else
	swapn8b(tp, *xpp, nelems);
# endif
	*xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_DOUBLE);
	return NC_NOERR;
}
#elif vax
int
ncmpix_getn_double_double(const void **xpp, MPI_Offset ndoubles, double *ip)
{
	double *const end = ip + ndoubles;

	while(ip < end)
	{
	struct vax_double *const vdp =
			 (struct vax_double *)ip;
	const struct ieee_double *const idp =
			 (const struct ieee_double *) (*xpp);
	{
		const struct dbl_limits *lim;
		int ii;
		for (ii = 0, lim = dbl_limits;
			ii < sizeof(dbl_limits)/sizeof(struct dbl_limits);
			ii++, lim++)
		{
			if ((idp->mant_lo == lim->ieee.mant_lo)
				&& (idp->mant_4 == lim->ieee.mant_4)
				&& (idp->mant_5 == lim->ieee.mant_5)
				&& (idp->mant_6 == lim->ieee.mant_6)
				&& (idp->exp_lo == lim->ieee.exp_lo)
				&& (idp->exp_hi == lim->ieee.exp_hi)
				)
			{
				*vdp = lim->d;
				goto doneit;
			}
		}
	}
	{
		unsigned exp = idp->exp_hi << 4 | idp->exp_lo;
		vdp->exp = exp - IEEE_DBL_BIAS + VAX_DBL_BIAS;
	}
	{
		unsigned mant_hi = ((idp->mant_6 << 16)
				 | (idp->mant_5 << 8)
				 | idp->mant_4);
		unsigned mant_lo = SWAP4(idp->mant_lo);
		vdp->mantissa1 = (mant_hi >> 13);
		vdp->mantissa2 = ((mant_hi & MASK(13)) << 3)
				| (mant_lo >> 29);
		vdp->mantissa3 = (mant_lo >> 13);
		vdp->mantissa4 = (mant_lo << 3);
	}
	doneit:
		vdp->sign = idp->sign;

		ip++;
		*xpp = (char *)(*xpp) + X_SIZEOF_DOUBLE;
	}
	return NC_NOERR;
}
	/* vax */

#else
int
ncmpix_getn_double_double(const void **xpp, MPI_Offset nelems, double *tp)
{
	const char *xp = *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		const int lstatus = ncmpix_get_double_double(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (const void *)xp;
	return status;
}

#endif

int
ncmpix_putn_double_schar(void **xpp, MPI_Offset nelems, const schar *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		int lstatus = ncmpix_put_double_schar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_double_uchar(void **xpp, MPI_Offset nelems, const uchar *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		int lstatus = ncmpix_put_double_uchar(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_double_short(void **xpp, MPI_Offset nelems, const short *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		int lstatus = ncmpix_put_double_short(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_double_int(void **xpp, MPI_Offset nelems, const int *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		int lstatus = ncmpix_put_double_int(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_double_long(void **xpp, MPI_Offset nelems, const long *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		int lstatus = ncmpix_put_double_long(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

int
ncmpix_putn_double_float(void **xpp, MPI_Offset nelems, const float *tp)
{
	char *xp = (char *) *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		int lstatus = ncmpix_put_double_float(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

#if X_SIZEOF_DOUBLE == SIZEOF_DOUBLE && !defined(NO_IEEE_FLOAT)
/* optimized version */
int
ncmpix_putn_double_double(void **xpp, MPI_Offset nelems, const double *tp)
{
# ifdef WORDS_BIGENDIAN
	(void) memcpy(*xpp, tp, nelems * X_SIZEOF_DOUBLE);
# else
	swapn8b(*xpp, tp, nelems);
# endif
	*xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_DOUBLE);
	return NC_NOERR;
}
#elif vax
int
ncmpix_putn_double_double(void **xpp, MPI_Offset ndoubles, const double *ip)
{
	const double *const end = ip + ndoubles;

	while(ip < end)
	{
	const struct vax_double *const vdp = 
			(const struct vax_double *)ip;
	struct ieee_double *const idp =
			 (struct ieee_double *) (*xpp);

	if ((vdp->mantissa4 > (dbl_limits[0].d.mantissa4 - 3)) &&
		(vdp->mantissa3 == dbl_limits[0].d.mantissa3) &&
		(vdp->mantissa2 == dbl_limits[0].d.mantissa2) &&
		(vdp->mantissa1 == dbl_limits[0].d.mantissa1) &&
		(vdp->exp == dbl_limits[0].d.exp))
	{
		*idp = dbl_limits[0].ieee;
		goto shipit;
	}
	if ((vdp->mantissa4 == dbl_limits[1].d.mantissa4) &&
		(vdp->mantissa3 == dbl_limits[1].d.mantissa3) &&
		(vdp->mantissa2 == dbl_limits[1].d.mantissa2) &&
		(vdp->mantissa1 == dbl_limits[1].d.mantissa1) &&
		(vdp->exp == dbl_limits[1].d.exp))
	{
		*idp = dbl_limits[1].ieee;
		goto shipit;
	}

	{
		unsigned exp = vdp->exp - VAX_DBL_BIAS + IEEE_DBL_BIAS;

		unsigned mant_lo = ((vdp->mantissa2 & MASK(3)) << 29) |
			(vdp->mantissa3 << 13) |
			((vdp->mantissa4 >> 3) & MASK(13));

		unsigned mant_hi = (vdp->mantissa1 << 13)
				 | (vdp->mantissa2 >> 3);

		if((vdp->mantissa4 & 7) > 4)
		{
			/* round up */
			mant_lo++;
			if(mant_lo == 0)
			{
				mant_hi++;
				if(mant_hi > 0xffffff)
				{
					mant_hi = 0;
					exp++;
				}
			}
		}

		idp->mant_lo = SWAP4(mant_lo);
		idp->mant_6 = mant_hi >> 16;
		idp->mant_5 = (mant_hi & 0xff00) >> 8;
		idp->mant_4 = mant_hi;
		idp->exp_hi = exp >> 4;
		idp->exp_lo = exp;
	}
		
	shipit:
		idp->sign = vdp->sign;

		ip++;
		*xpp = (char *)(*xpp) + X_SIZEOF_DOUBLE;
	}
	return NC_NOERR;
}
	/* vax */

#else
int
ncmpix_putn_double_double(void **xpp, MPI_Offset nelems, const double *tp)
{
	char *xp = *xpp;
	int status = NC_NOERR;

	for( ; nelems != 0; nelems--, xp += X_SIZEOF_DOUBLE, tp++)
	{
		int lstatus = ncmpix_put_double_double(xp, tp);
		if(lstatus != NC_NOERR)
			status = lstatus;
	}

	*xpp = (void *)xp;
	return status;
}

#endif


/*
 * Other aggregate conversion functions.
 */

/* text */

int
ncmpix_getn_text(const void **xpp, MPI_Offset nelems, char *tp)
{
	(void) memcpy(tp, *xpp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);
	return NC_NOERR;

}

int
ncmpix_pad_getn_text(const void **xpp, MPI_Offset nelems, char *tp)
{
	MPI_Offset rndup = nelems % X_ALIGN;

	if(rndup)
		rndup = X_ALIGN - rndup;

	(void) memcpy(tp, *xpp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems + rndup);

	return NC_NOERR;

}

int
ncmpix_putn_text(void **xpp, MPI_Offset nelems, const char *tp)
{
	(void) memcpy(*xpp, tp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);

	return NC_NOERR;

}

int
ncmpix_pad_putn_text(void **xpp, MPI_Offset nelems, const char *tp)
{
	MPI_Offset rndup = nelems % X_ALIGN;

	if(rndup)
		rndup = X_ALIGN - rndup;

	(void) memcpy(*xpp, tp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);

	if(rndup)
	{
		(void) memcpy(*xpp, nada, rndup);
		*xpp = (void *)((char *)(*xpp) + rndup);
	}
	
	return NC_NOERR;

}


/* opaque */

int
ncmpix_getn_void(const void **xpp, MPI_Offset nelems, void *tp)
{
	(void) memcpy(tp, *xpp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);
	return NC_NOERR;

}

int
ncmpix_pad_getn_void(const void **xpp, MPI_Offset nelems, void *tp)
{
	MPI_Offset rndup = nelems % X_ALIGN;

	if(rndup)
		rndup = X_ALIGN - rndup;

	(void) memcpy(tp, *xpp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems + rndup);

	return NC_NOERR;

}

int
ncmpix_putn_void(void **xpp, MPI_Offset nelems, const void *tp)
{
	(void) memcpy(*xpp, tp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);

	return NC_NOERR;

}

int
ncmpix_pad_putn_void(void **xpp, MPI_Offset nelems, const void *tp)
{
	MPI_Offset rndup = nelems % X_ALIGN;

	if(rndup)
		rndup = X_ALIGN - rndup;

	(void) memcpy(*xpp, tp, nelems);
	*xpp = (void *)((char *)(*xpp) + nelems);

	if(rndup)
	{
		(void) memcpy(*xpp, nada, rndup);
		*xpp = (void *)((char *)(*xpp) + rndup);
	}
	
	return NC_NOERR;

}

