/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: rnd.h 816 2010-05-14 23:34:49Z wkliao $ */
#ifndef _RNDUP

/* useful for aligning memory */
#define	_RNDUP(x, unit)  ((((x) + (unit) - 1) / (unit)) \
	* (unit))
#define	_RNDDOWN(x, unit)  ((x) - ((x)%(unit)))

#define M_RND_UNIT	(sizeof(double))
#define	M_RNDUP(x) _RNDUP(x, M_RND_UNIT)
#define	M_RNDDOWN(x)  __RNDDOWN(x, M_RND_UNIT)

#endif
