/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: attr.c 829 2010-05-26 20:17:57Z wkliao $ */

#include "nc.h"
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h>
#include <assert.h>
#include "nc.h"
#include "ncx.h"
#include "fbits.h"
#include "rnd.h"
#include "macro.h"


/*
 * Free attr
 * Formerly
NC_free_attr()
 */
void
ncmpii_free_NC_attr(NC_attr *attrp)
{

	if(attrp == NULL)
		return;
	ncmpii_free_NC_string(attrp->name);
	NCI_Free(attrp);
}


/*
 * How much space will 'nelems' of 'type' take in
 *  external representation (as the values of an attribute)?
 */
static size_t
ncmpix_len_NC_attrV(nc_type type, size_t nelems)
{
	switch(type) {
	case NC_BYTE:
	case NC_CHAR:
		return ncmpix_len_char(nelems);
	case NC_SHORT:
		return ncmpix_len_short(nelems);
	case NC_INT:
		return ncmpix_len_int(nelems);
	case NC_FLOAT:
		return ncmpix_len_float(nelems);
	case NC_DOUBLE:
		return ncmpix_len_double(nelems);
	default:
		assert("ncmpix_len_NC_attr bad type" == 0);
	}
	return 0;
}


NC_attr *
ncmpii_new_x_NC_attr(
	NC_string *strp,
	nc_type type,
	MPI_Offset nelems)
{
	NC_attr *attrp;
	const size_t xsz = ncmpix_len_NC_attrV(type, nelems);
	size_t sz = M_RNDUP(sizeof(NC_attr));

	assert(!(xsz == 0 && nelems != 0));

	sz += xsz;

	attrp = (NC_attr *) NCI_Malloc(sz);
	if(attrp == NULL )
		return NULL;

	attrp->xsz = xsz;

	attrp->name = strp;
	attrp->type = type;
	attrp->nelems = nelems;
	if(xsz != 0)
		attrp->xvalue = (char *)attrp + M_RNDUP(sizeof(NC_attr));
	else
		attrp->xvalue = NULL;

	return(attrp);
}


/*
 * Formerly
NC_new_attr(name,type,count,value)
 */
static NC_attr *
ncmpii_new_NC_attr(
	const char *name,
	nc_type type,
	MPI_Offset nelems)
{
	NC_string *strp;
	NC_attr *attrp;

	assert(name != NULL && *name != 0);

	strp = ncmpii_new_NC_string(strlen(name), name);
	if(strp == NULL)
		return NULL;
	
	attrp = ncmpii_new_x_NC_attr(strp, type, nelems);
	if(attrp == NULL)
	{
		ncmpii_free_NC_string(strp);
		return NULL;
	}

	return(attrp);
}


static NC_attr *
dup_NC_attr(const NC_attr *rattrp)
{
	NC_attr *attrp = ncmpii_new_NC_attr(rattrp->name->cp,
		 rattrp->type, rattrp->nelems);
	if(attrp == NULL)
		return NULL;
	(void) memcpy(attrp->xvalue, rattrp->xvalue, rattrp->xsz);
	return attrp;
}

/* attrarray */

/*
 * Free the stuff "in" (referred to by) an NC_attrarray.
 * Leaves the array itself allocated.
 */
void
ncmpii_free_NC_attrarrayV0(NC_attrarray *ncap)
{
	assert(ncap != NULL);

	if(ncap->nelems == 0)
		return;

	assert(ncap->value != NULL);

	{
		NC_attr **app = ncap->value;
		NC_attr *const *const end = &app[ncap->nelems];
		for( /*NADA*/; app < end; app++)
		{
			ncmpii_free_NC_attr(*app);
			*app = NULL;
		}
	}
	ncap->nelems = 0;
}


/*
 * Free NC_attrarray values.
 * formerly
NC_free_array()
 */
void
ncmpii_free_NC_attrarrayV(NC_attrarray *ncap)
{
	assert(ncap != NULL);
	
	if(ncap->nalloc == 0)
		return;

	assert(ncap->value != NULL);

	ncmpii_free_NC_attrarrayV0(ncap);

	NCI_Free(ncap->value);
	ncap->value = NULL;
	ncap->nalloc = 0;
}


int
ncmpii_dup_NC_attrarrayV(NC_attrarray *ncap, const NC_attrarray *ref)
{
	int status = NC_NOERR;

	assert(ref != NULL);
	assert(ncap != NULL);

	if(ref->nelems != 0)
	{
		const size_t sz = ref->nelems * sizeof(NC_attr *);
		ncap->value = (NC_attr **) NCI_Malloc(sz);
		if(ncap->value == NULL)
			return NC_ENOMEM;

		(void) memset(ncap->value, 0, sz);
		ncap->nalloc = ref->nelems;
	}

	ncap->nelems = 0;
	{
		NC_attr **app = ncap->value;
		const NC_attr **drpp = (const NC_attr **)ref->value;
		NC_attr *const *const end = &app[ref->nelems];
		for( /*NADA*/; app < end; drpp++, app++, ncap->nelems++)
		{
			*app = dup_NC_attr(*drpp);
			if(*app == NULL)
			{
				status = NC_ENOMEM;
				break;
			}
		}
	}

	if(status != NC_NOERR)
	{
		ncmpii_free_NC_attrarrayV(ncap);
		return status;
	}

	assert(ncap->nelems == ref->nelems);

	return NC_NOERR;
}


/*
 * Add a new handle on the end of an array of handles
 * Formerly
NC_incr_array(array, tail)
 */
static int
incr_NC_attrarray(NC_attrarray *ncap, NC_attr *newelemp)
{
	NC_attr **vp;

	assert(ncap != NULL);

	if(ncap->nalloc == 0)
	{
		assert(ncap->nelems == 0);
		vp = (NC_attr **) NCI_Malloc(NC_ARRAY_GROWBY * sizeof(NC_attr *));
		if(vp == NULL)
			return NC_ENOMEM;

		ncap->value = vp;
		ncap->nalloc = NC_ARRAY_GROWBY;
	}
	else if(ncap->nelems +1 > ncap->nalloc)
	{
		vp = (NC_attr **) NCI_Realloc(ncap->value,
			(ncap->nalloc + NC_ARRAY_GROWBY) * sizeof(NC_attr *));
		if(vp == NULL)
			return NC_ENOMEM;
	
		ncap->value = vp;
		ncap->nalloc += NC_ARRAY_GROWBY;
	}

	if(newelemp != NULL)
	{
		ncap->value[ncap->nelems] = newelemp;
		ncap->nelems++;
	}
	return NC_NOERR;
}


static NC_attr *
elem_NC_attrarray(const NC_attrarray *ncap, MPI_Offset elem)
{
	assert(ncap != NULL);
	if((elem < 0) || ncap->nelems == 0 || elem >= ncap->nelems)
		return NULL;

	assert(ncap->value != NULL);

	return ncap->value[elem];
}

/* End attarray per se */

/*
 * Given ncp and varid, return ptr to array of attributes
 *  else NULL on error
 */
static NC_attrarray *
NC_attrarray0( NC *ncp, int varid)
{
	NC_attrarray *ap;

	if(varid == NC_GLOBAL) /* Global attribute, attach to cdf */
	{
		ap = &ncp->attrs;
	}
	else if(varid >= 0 && (size_t) varid < ncp->vars.nelems)
	{
		NC_var **vpp;
		vpp = (NC_var **)ncp->vars.value;
		vpp += varid;
		ap = &(*vpp)->attrs;
	} else {
		ap = NULL;
	}
	return(ap);
}


/*
 * Step thru NC_ATTRIBUTE array, seeking match on name.
 *  return match or NULL if Not Found.
 */
NC_attr **
ncmpii_NC_findattr(const NC_attrarray *ncap, const char *name)
{
	NC_attr **attrpp;
	size_t attrid;
	size_t slen;

	assert(ncap != NULL);

	if(ncap->nelems == 0)
		return NULL;

	attrpp = (NC_attr **) ncap->value;

	slen = strlen(name);

	for(attrid = 0; attrid < ncap->nelems; attrid++, attrpp++)
	{
		if(strlen((*attrpp)->name->cp) == slen &&
			strncmp((*attrpp)->name->cp, name, slen) == 0)
		{
			return(attrpp); /* Normal return */
		}
	}
	return(NULL);
}


/*
 * Look up by ncid, varid and name, return NULL if not found
 */
static int 
NC_lookupattr(int ncid,
	int varid,
	const char *name, /* attribute name */
	NC_attr **attrpp) /* modified on return */
{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **tmp;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	tmp = ncmpii_NC_findattr(ncap, name);
	if(tmp == NULL)
		return NC_ENOTATT;

	if(attrpp != NULL)
		*attrpp = *tmp;

	return NC_NOERR;
}

/* Public */

int
ncmpi_inq_attname(int ncid, int varid, int attnum, char *name)

{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr *attrp;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	attrp = elem_NC_attrarray(ncap, attnum);
	if(attrp == NULL)
		return NC_ENOTATT;

	(void) strncpy(name, attrp->name->cp, attrp->name->nchars);
	name[attrp->name->nchars] = 0;

	return NC_NOERR;
}


int 
ncmpi_inq_attid(int ncid, int varid, const char *name, int *attnump)
{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;
	

	attrpp = ncmpii_NC_findattr(ncap, name);
	if(attrpp == NULL)
		return NC_ENOTATT;

	if(attnump != NULL)
		*attnump = (int)(attrpp - ncap->value);

	return NC_NOERR;
}

int 
ncmpi_inq_atttype(int ncid, int varid, const char *name, nc_type *datatypep)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(datatypep != NULL)
		*datatypep = attrp->type;

	return NC_NOERR;
}

int 
ncmpi_inq_attlen(int ncid, int varid, const char *name, MPI_Offset *lenp)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(lenp != NULL)
		*lenp = attrp->nelems;

	return NC_NOERR;
}

int
ncmpi_inq_att(int ncid,
	int varid,
	const char *name, /* input, attribute name */
	nc_type *datatypep,
	MPI_Offset *lenp)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(datatypep != NULL)
		*datatypep = attrp->type;
	if(lenp != NULL)
		*lenp = attrp->nelems;

	return NC_NOERR;
}


int
ncmpi_rename_att( int ncid, int varid, const char *name, const char *newname)
{
    int status;
    NC *ncp;
    NC_attrarray *ncap;
    NC_attr **tmp;
    NC_attr *attrp;
    NC_string *newStr, *old;

    /* sortof inline clone of NC_lookupattr() */
    status = ncmpii_NC_check_id(ncid, &ncp);
    if (status != NC_NOERR)
        return status;

    if (NC_readonly(ncp))
        return NC_EPERM;

    ncap = NC_attrarray0(ncp, varid);
    if (ncap == NULL)
        return NC_ENOTVAR;

/* bugs found by Jianwei Li
        status = ncmpii_NC_check_name(name);
*/
    status = ncmpii_NC_check_name(newname);
    if (status != NC_NOERR)
        return status;

    tmp = ncmpii_NC_findattr(ncap, name);
    if (tmp == NULL)
        return NC_ENOTATT;
    attrp = *tmp;
    /* end inline clone NC_lookupattr() */

    if (ncmpii_NC_findattr(ncap, newname) != NULL) {
        /* name in use */
        return NC_ENAMEINUSE;
    }

    old = attrp->name;
    if (NC_indef(ncp)) {
        newStr = ncmpii_new_NC_string(strlen(newname), newname);
        if (newStr == NULL)
            return NC_ENOMEM;
        attrp->name = newStr;
        ncmpii_free_NC_string(old);
        return NC_NOERR;
    }
    /* else, no in define mode */
    status = ncmpii_set_NC_string(old, newname);
    if (status != NC_NOERR)
        return status;

    set_NC_hdirty(ncp);

    if (NC_doHsync(ncp)) {
        status = ncmpii_NC_sync(ncp, 1);
        if (status != NC_NOERR)
            return status;
    }

    return NC_NOERR;
}


int
ncmpi_copy_att(int ncid_in, int varid_in, const char *name, int ncid_out, int ovarid)
{
	int status;
	NC_attr *iattrp;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;
	NC_attr *old = NULL;
	NC_attr *attrp;

	status = NC_lookupattr(ncid_in, varid_in, name, &iattrp);
	if(status != NC_NOERR)
		return status;

	status = ncmpii_NC_check_id(ncid_out, &ncp);
	if(status != NC_NOERR)
		return status;

	if(NC_readonly(ncp))
		return NC_EPERM;

	ncap = NC_attrarray0(ncp, ovarid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	attrpp = ncmpii_NC_findattr(ncap, name);
	if(attrpp != NULL) /* name in use */
	{
		if(!NC_indef(ncp) )
		{
			attrp = *attrpp; /* convenience */
	
			if(iattrp->xsz > attrp->xsz)
				return NC_ENOTINDEFINE;
			/* else, we can reuse existing without redef */
			
			attrp->xsz = iattrp->xsz;
			attrp->type = iattrp->type;
			attrp->nelems = iattrp->nelems;

			(void) memcpy(attrp->xvalue, iattrp->xvalue,
				iattrp->xsz);
			
			set_NC_hdirty(ncp);

			if(NC_doHsync(ncp))
			{
				status = ncmpii_NC_sync(ncp, 1);
				if(status != NC_NOERR)
					return status;
			}

			return NC_NOERR;
		}
		/* else, redefine using existing array slot */
		old = *attrpp;
	} 
	else
	{
		if(!NC_indef(ncp))
			return NC_ENOTINDEFINE;

		if(ncap->nelems >= NC_MAX_ATTRS)
			return NC_EMAXATTS;
	}

	attrp = ncmpii_new_NC_attr(name, iattrp->type, iattrp->nelems);
	if(attrp == NULL)
		return NC_ENOMEM;

	(void) memcpy(attrp->xvalue, iattrp->xvalue,
		iattrp->xsz);

	if(attrpp != NULL)
	{
		assert(old != NULL);
		*attrpp = attrp;
		ncmpii_free_NC_attr(old);
	}
	else
	{
		status = incr_NC_attrarray(ncap, attrp);
		if(status != NC_NOERR)
		{
			ncmpii_free_NC_attr(attrp);
			return status;
		}
	}

	return NC_NOERR;
}


int
ncmpi_del_att(int ncid, int varid, const char *name)
{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;
	NC_attr *old = NULL;
	int attrid;
	MPI_Offset slen;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	if(!NC_indef(ncp))
		return NC_ENOTINDEFINE;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

			/* sortof inline ncmpii_NC_findattr() */
	slen = strlen(name);

	attrpp = (NC_attr **) ncap->value;
	for(attrid = 0; (size_t) attrid < ncap->nelems; attrid++, attrpp++)
	{
		if( slen == (*attrpp)->name->nchars &&
			strncmp(name, (*attrpp)->name->cp, slen) == 0)
		{
			old = *attrpp;
			break;
		}
	}
	if( (size_t) attrid == ncap->nelems )
		return NC_ENOTATT;
			/* end inline NC_findattr() */

	/* shuffle down */
	for(attrid++; (size_t) attrid < ncap->nelems; attrid++)
	{
		*attrpp = *(attrpp + 1);
		attrpp++;
	}
	*attrpp = NULL;
	/* decrement count */
	ncap->nelems--;

	ncmpii_free_NC_attr(old);

	return NC_NOERR;
}


static int
ncmpix_pad_putn_Iuchar(void **xpp, MPI_Offset nelems, const uchar *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_putn_schar_uchar(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_putn_short_uchar(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_putn_int_uchar(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_putn_float_uchar(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_putn_double_uchar(xpp, nelems, tp);
	default: 
		assert("ncmpix_pad_putn_Iuchar invalid type" == 0);
		return NC_EBADTYPE;
	}
}

static int
ncmpix_pad_getn_Iuchar(const void **xpp, MPI_Offset nelems, uchar *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_getn_schar_uchar(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_getn_short_uchar(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_getn_int_uchar(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_getn_float_uchar(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_getn_double_uchar(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_getn_Iuchar invalid type" == 0);
		return NC_EBADTYPE;
	}
}


static int
ncmpix_pad_putn_Ischar(void **xpp, MPI_Offset nelems, const schar *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_putn_schar_schar(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_putn_short_schar(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_putn_int_schar(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_putn_float_schar(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_putn_double_schar(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_putn_Ischar invalid type" == 0);
		return NC_EBADTYPE;
	}
}

static int
ncmpix_pad_getn_Ischar(const void **xpp, MPI_Offset nelems, schar *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_getn_schar_schar(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_getn_short_schar(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_getn_int_schar(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_getn_float_schar(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_getn_double_schar(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_getn_Ischar invalid type" == 0);
		return NC_EBADTYPE;
	}
}


static int
ncmpix_pad_putn_Ishort(void **xpp, MPI_Offset nelems, const short *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_putn_schar_short(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_putn_short_short(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_putn_int_short(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_putn_float_short(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_putn_double_short(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_putn_Ishort invalid type" == 0);
		return NC_EBADTYPE;
	}
}

static int
ncmpix_pad_getn_Ishort(const void **xpp, MPI_Offset nelems, short *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_getn_schar_short(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_getn_short_short(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_getn_int_short(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_getn_float_short(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_getn_double_short(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_getn_Ishort invalid type" == 0);
		return NC_EBADTYPE;
	}
}


static int
ncmpix_pad_putn_Iint(void **xpp, MPI_Offset nelems, const int *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_putn_schar_int(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_putn_short_int(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_putn_int_int(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_putn_float_int(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_putn_double_int(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_putn_Iint invalid type" == 0);
		return NC_EBADTYPE;
	}
}

static int
ncmpix_pad_getn_Iint(const void **xpp, MPI_Offset nelems, int *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_getn_schar_int(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_getn_short_int(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_getn_int_int(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_getn_float_int(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_getn_double_int(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_getn_Iint invalid type" == 0);
		return NC_EBADTYPE;
	}
}


static int
ncmpix_pad_putn_Ilong(void **xpp, MPI_Offset nelems, const long *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_putn_schar_long(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_putn_short_long(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_putn_int_long(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_putn_float_long(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_putn_double_long(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_putn_Ilong invalid type" == 0);
		return NC_EBADTYPE;
	}
}

static int
ncmpix_pad_getn_Ilong(const void **xpp, MPI_Offset nelems, long *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_getn_schar_long(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_getn_short_long(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_getn_int_long(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_getn_float_long(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_getn_double_long(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_getn_Ilong invalid type" == 0);
		return NC_EBADTYPE;
	}
}


static int
ncmpix_pad_putn_Ifloat(void **xpp, MPI_Offset nelems, const float *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_putn_schar_float(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_putn_short_float(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_putn_int_float(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_putn_float_float(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_putn_double_float(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_putn_Ifloat invalid type" == 0);
		return NC_EBADTYPE;
	}
}

static int
ncmpix_pad_getn_Ifloat(const void **xpp, MPI_Offset nelems, float *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_getn_schar_float(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_getn_short_float(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_getn_int_float(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_getn_float_float(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_getn_double_float(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_getn_Ifloat invalid type" == 0);
		return NC_EBADTYPE;
	}
}


static int
ncmpix_pad_putn_Idouble(void **xpp, MPI_Offset nelems, const double *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_putn_schar_double(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_putn_short_double(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_putn_int_double(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_putn_float_double(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_putn_double_double(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_putn_Idouble invalid type" == 0);
		return NC_EBADTYPE;
	}
}

static int
ncmpix_pad_getn_Idouble(const void **xpp, MPI_Offset nelems, double *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return ncmpix_pad_getn_schar_double(xpp, nelems, tp);
	case NC_SHORT:
		return ncmpix_pad_getn_short_double(xpp, nelems, tp);
	case NC_INT:
		return ncmpix_getn_int_double(xpp, nelems, tp);
	case NC_FLOAT:
		return ncmpix_getn_float_double(xpp, nelems, tp);
	case NC_DOUBLE:
		return ncmpix_getn_double_double(xpp, nelems, tp);
	default:
		assert("ncmpix_pad_getn_Idouble invalid type" == 0);
		return NC_EBADTYPE;
	}
}



int
ncmpi_put_att_text(int ncid, int varid, const char *name,
	MPI_Offset nelems, const char *value)
{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;
	NC_attr *old = NULL;
	NC_attr *attrp;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	if(NC_readonly(ncp))
		return NC_EPERM;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	status = ncmpii_NC_check_name(name);
	if(status != NC_NOERR)
		return status;

	if(nelems < 0 || nelems > X_INT_MAX) /* backward compat */
		return NC_EINVAL; /* Invalid nelems */

	if(nelems != 0 && value == NULL)
		return NC_EINVAL; /* Null arg */

	attrpp = ncmpii_NC_findattr(ncap, name);
	if(attrpp != NULL) /* name in use */
	{
		if(!NC_indef(ncp) )
		{
			const size_t xsz = ncmpix_len_NC_attrV(NC_CHAR, nelems);
			attrp = *attrpp; /* convenience */
	
			if(xsz > attrp->xsz)
				return NC_ENOTINDEFINE;
			/* else, we can reuse existing without redef */
			
			attrp->xsz = xsz;
			attrp->type = NC_CHAR;
			attrp->nelems = nelems;

			if(nelems != 0)
			{
				void *xp = attrp->xvalue;
				status = ncmpix_pad_putn_text(&xp, nelems, value);
				if(status != NC_NOERR)
					return status;
			}
			
			set_NC_hdirty(ncp);

			if(NC_doHsync(ncp))
			{
				status = ncmpii_NC_sync(ncp, 1);
				if(status != NC_NOERR)
					return status;
			}

			return NC_NOERR;
		}
		/* else, redefine using existing array slot */
		old = *attrpp;
	} 
	else
	{
		if(!NC_indef(ncp))
			return NC_ENOTINDEFINE;

		if(ncap->nelems >= NC_MAX_ATTRS)
			return NC_EMAXATTS;
	}

	attrp = ncmpii_new_NC_attr(name, NC_CHAR, nelems);
	if(attrp == NULL)
		return NC_ENOMEM;

	if(nelems != 0)
	{
		void *xp = attrp->xvalue;
		status = ncmpix_pad_putn_text(&xp, nelems, value);
		if(status != NC_NOERR)
			return status;
	}

	if(attrpp != NULL)
	{
		assert(old != NULL);
		*attrpp = attrp;
		ncmpii_free_NC_attr(old);
	}
	else
	{
		status = incr_NC_attrarray(ncap, attrp);
		if(status != NC_NOERR)
		{
			ncmpii_free_NC_attr(attrp);
			return status;
		}
	}

	return NC_NOERR;
}


int
ncmpi_get_att_text(int ncid, int varid, const char *name, char *str)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(attrp->nelems == 0)
		return NC_NOERR;

	if(attrp->type != NC_CHAR)
		return NC_ECHAR;

	/* else */
	{
		const void *xp = attrp->xvalue;
		return ncmpix_pad_getn_text(&xp, attrp->nelems, str);
	}
}




int
ncmpi_put_att_schar(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const signed char *value)
{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;
	NC_attr *old = NULL;
	NC_attr *attrp;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	if(NC_readonly(ncp))
		return NC_EPERM;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	status = ncmpii_cktype(type);
	if(status != NC_NOERR)
		return status;

	if(type == NC_CHAR)
		return NC_ECHAR;

	if(nelems < 0 || nelems > X_INT_MAX) /* backward compat */
		return NC_EINVAL; /* Invalid nelems */

	if(nelems != 0 && value == NULL)
		return NC_EINVAL; /* Null arg */

	attrpp = ncmpii_NC_findattr(ncap, name);
	if(attrpp != NULL) /* name in use */
	{
		if(!NC_indef(ncp) )
		{
			const size_t xsz = ncmpix_len_NC_attrV(type, nelems);
			attrp = *attrpp; /* convenience */
	
			if(xsz > attrp->xsz)
				return NC_ENOTINDEFINE;
			/* else, we can reuse existing without redef */
			
			attrp->xsz = xsz;
			attrp->type = type;
			attrp->nelems = nelems;

			if(nelems != 0)
			{
				void *xp = attrp->xvalue;
				status = ncmpix_pad_putn_Ischar(&xp, nelems,
					value, type);
			}
			
			set_NC_hdirty(ncp);

			if(NC_doHsync(ncp))
			{
				const int lstatus = ncmpii_NC_sync(ncp, 1);
				/*
				 * N.B.: potentially overrides NC_ERANGE
				 * set by ncmpix_pad_putn_Ischar
				 */
				if(lstatus != NC_NOERR)
					return lstatus;
			}

			return status;
		}
		/* else, redefine using existing array slot */
		old = *attrpp;
	} 
	else
	{
		if(!NC_indef(ncp))
			return NC_ENOTINDEFINE;

		if(ncap->nelems >= NC_MAX_ATTRS)
			return NC_EMAXATTS;
	}

	status = ncmpii_NC_check_name(name);
	if(status != NC_NOERR)
		return status;

	attrp = ncmpii_new_NC_attr(name, type, nelems);
	if(attrp == NULL)
		return NC_ENOMEM;

	if(nelems != 0)
	{
		void *xp = attrp->xvalue;
		status = ncmpix_pad_putn_Ischar(&xp, nelems,
			value, type);
	}

	if(attrpp != NULL)
	{
		assert(old != NULL);
		*attrpp = attrp;
		ncmpii_free_NC_attr(old);
	}
	else
	{
		const int lstatus = incr_NC_attrarray(ncap, attrp);
		/*
		 * N.B.: potentially overrides NC_ERANGE
		 * set by ncmpix_pad_putn_Ischar
		 */
		if(lstatus != NC_NOERR)
		{
			ncmpii_free_NC_attr(attrp);
			return lstatus;
		}
	}

	return status;
}

int
ncmpi_get_att_schar(int ncid, int varid, const char *name, signed char *tp)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(attrp->nelems == 0)
		return NC_NOERR;

	if(attrp->type == NC_CHAR)
		return NC_ECHAR;

	{
	const void *xp = attrp->xvalue;
	return ncmpix_pad_getn_Ischar(&xp, attrp->nelems, tp, attrp->type);
	}
}


int
ncmpi_put_att_uchar(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const unsigned char *value)
{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;
	NC_attr *old = NULL;
	NC_attr *attrp;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	if(NC_readonly(ncp))
		return NC_EPERM;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	status = ncmpii_cktype(type);
	if(status != NC_NOERR)
		return status;

	if(type == NC_CHAR)
		return NC_ECHAR;

	if( (nelems < 0) || nelems > X_INT_MAX) /* backward compat */
		return NC_EINVAL; /* Invalid nelems */

	if(nelems != 0 && value == NULL)
		return NC_EINVAL; /* Null arg */

	attrpp = ncmpii_NC_findattr(ncap, name);
	if(attrpp != NULL) /* name in use */
	{
		if(!NC_indef(ncp) )
		{
			const size_t xsz = ncmpix_len_NC_attrV(type, nelems);
			attrp = *attrpp; /* convenience */
	
			if(xsz > attrp->xsz)
				return NC_ENOTINDEFINE;
			/* else, we can reuse existing without redef */
			
			attrp->xsz = xsz;
			attrp->type = type;
			attrp->nelems = nelems;

			if(nelems != 0)
			{
				void *xp = attrp->xvalue;
				status = ncmpix_pad_putn_Iuchar(&xp, nelems,
					value, type);
			}
			
			set_NC_hdirty(ncp);

			if(NC_doHsync(ncp))
			{
				const int lstatus = ncmpii_NC_sync(ncp, 1);
				/*
				 * N.B.: potentially overrides NC_ERANGE
				 * set by ncmpix_pad_putn_Iuchar
				 */
				if(lstatus != NC_NOERR)
					return lstatus;
			}

			return status;
		}
		/* else, redefine using existing array slot */
		old = *attrpp;
	} 
	else
	{
		if(!NC_indef(ncp))
			return NC_ENOTINDEFINE;

		if(ncap->nelems >= NC_MAX_ATTRS)
			return NC_EMAXATTS;
	}

	status = ncmpii_NC_check_name(name);
	if(status != NC_NOERR)
		return status;

	attrp = ncmpii_new_NC_attr(name, type, nelems);
	if(attrp == NULL)
		return NC_ENOMEM;

	if(nelems != 0)
	{
		void *xp = attrp->xvalue;
		status = ncmpix_pad_putn_Iuchar(&xp, nelems,
			value, type);
	}

	if(attrpp != NULL)
	{
		assert(old != NULL);
		*attrpp = attrp;
		ncmpii_free_NC_attr(old);
	}
	else
	{
		const int lstatus = incr_NC_attrarray(ncap, attrp);
		/*
		 * N.B.: potentially overrides NC_ERANGE
		 * set by ncmpix_pad_putn_Iuchar
		 */
		if(lstatus != NC_NOERR)
		{
			ncmpii_free_NC_attr(attrp);
			return lstatus;
		}
	}

	return status;
}

int
ncmpi_get_att_uchar(int ncid, int varid, const char *name, unsigned char *tp)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(attrp->nelems == 0)
		return NC_NOERR;

	if(attrp->type == NC_CHAR)
		return NC_ECHAR;

	{
	const void *xp = attrp->xvalue;
	return ncmpix_pad_getn_Iuchar(&xp, attrp->nelems, tp, attrp->type);
	}
}


int
ncmpi_put_att_short(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const short *value)
{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;
	NC_attr *old = NULL;
	NC_attr *attrp;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	if(NC_readonly(ncp))
		return NC_EPERM;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	status = ncmpii_cktype(type);
	if(status != NC_NOERR)
		return status;

	if(type == NC_CHAR)
		return NC_ECHAR;

	if( (nelems < 0) || (nelems > X_INT_MAX)) /* backward compat */
		return NC_EINVAL; /* Invalid nelems */

	if(nelems != 0 && value == NULL)
		return NC_EINVAL; /* Null arg */

	attrpp = ncmpii_NC_findattr(ncap, name);
	if(attrpp != NULL) /* name in use */
	{
		if(!NC_indef(ncp) )
		{
			const size_t xsz = ncmpix_len_NC_attrV(type, nelems);
			attrp = *attrpp; /* convenience */
	
			if(xsz > attrp->xsz)
				return NC_ENOTINDEFINE;
			/* else, we can reuse existing without redef */
			
			attrp->xsz = xsz;
			attrp->type = type;
			attrp->nelems = nelems;

			if(nelems != 0)
			{
				void *xp = attrp->xvalue;
				status = ncmpix_pad_putn_Ishort(&xp, nelems,
					value, type);
			}
			
			set_NC_hdirty(ncp);

			if(NC_doHsync(ncp))
			{
				const int lstatus = ncmpii_NC_sync(ncp, 1);
				/*
				 * N.B.: potentially overrides NC_ERANGE
				 * set by ncmpix_pad_putn_Ishort
				 */
				if(lstatus != NC_NOERR)
					return lstatus;
			}

			return status;
		}
		/* else, redefine using existing array slot */
		old = *attrpp;
	} 
	else
	{
		if(!NC_indef(ncp))
			return NC_ENOTINDEFINE;

		if(ncap->nelems >= NC_MAX_ATTRS)
			return NC_EMAXATTS;
	}

	status = ncmpii_NC_check_name(name);
	if(status != NC_NOERR)
		return status;

	attrp = ncmpii_new_NC_attr(name, type, nelems);
	if(attrp == NULL)
		return NC_ENOMEM;

	if(nelems != 0)
	{
		void *xp = attrp->xvalue;
		status = ncmpix_pad_putn_Ishort(&xp, nelems,
			value, type);
	}

	if(attrpp != NULL)
	{
		assert(old != NULL);
		*attrpp = attrp;
		ncmpii_free_NC_attr(old);
	}
	else
	{
		const int lstatus = incr_NC_attrarray(ncap, attrp);
		/*
		 * N.B.: potentially overrides NC_ERANGE
		 * set by ncmpix_pad_putn_Ishort
		 */
		if(lstatus != NC_NOERR)
		{
			ncmpii_free_NC_attr(attrp);
			return lstatus;
		}
	}

	return status;
}

int
ncmpi_get_att_short(int ncid, int varid, const char *name, short *tp)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(attrp->nelems == 0)
		return NC_NOERR;

	if(attrp->type == NC_CHAR)
		return NC_ECHAR;

	{
	const void *xp = attrp->xvalue;
	return ncmpix_pad_getn_Ishort(&xp, attrp->nelems, tp, attrp->type);
	}
}


int
ncmpi_put_att_int(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const int *value)
{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;
	NC_attr *old = NULL;
	NC_attr *attrp;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	if(NC_readonly(ncp))
		return NC_EPERM;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	status = ncmpii_cktype(type);
	if(status != NC_NOERR)
		return status;

	if(type == NC_CHAR)
		return NC_ECHAR;

	if((nelems < 0) ||  (nelems > X_INT_MAX) )/* backward compat */
		return NC_EINVAL; /* Invalid nelems */

	if(nelems != 0 && value == NULL)
		return NC_EINVAL; /* Null arg */

	attrpp = ncmpii_NC_findattr(ncap, name);
	if(attrpp != NULL) /* name in use */
	{
		if(!NC_indef(ncp) )
		{
			const size_t xsz = ncmpix_len_NC_attrV(type, nelems);
			attrp = *attrpp; /* convenience */
	
			if(xsz > attrp->xsz)
				return NC_ENOTINDEFINE;
			/* else, we can reuse existing without redef */
			
			attrp->xsz = xsz;
			attrp->type = type;
			attrp->nelems = nelems;

			if(nelems != 0)
			{
				void *xp = attrp->xvalue;
				status = ncmpix_pad_putn_Iint(&xp, nelems,
					value, type);
			}
			
			set_NC_hdirty(ncp);

			if(NC_doHsync(ncp))
			{
				const int lstatus = ncmpii_NC_sync(ncp, 1);
				/*
				 * N.B.: potentially overrides NC_ERANGE
				 * set by ncmpix_pad_putn_Iint
				 */
				if(lstatus != NC_NOERR)
					return lstatus;
			}

			return status;
		}
		/* else, redefine using existing array slot */
		old = *attrpp;
	} 
	else
	{
		if(!NC_indef(ncp))
			return NC_ENOTINDEFINE;

		if(ncap->nelems >= NC_MAX_ATTRS)
			return NC_EMAXATTS;
	}

	status = ncmpii_NC_check_name(name);
	if(status != NC_NOERR)
		return status;

	attrp = ncmpii_new_NC_attr(name, type, nelems);
	if(attrp == NULL)
		return NC_ENOMEM;

	if(nelems != 0)
	{
		void *xp = attrp->xvalue;
		status = ncmpix_pad_putn_Iint(&xp, nelems,
			value, type);
	}

	if(attrpp != NULL)
	{
		assert(old != NULL);
		*attrpp = attrp;
		ncmpii_free_NC_attr(old);
	}
	else
	{
		const int lstatus = incr_NC_attrarray(ncap, attrp);
		/*
		 * N.B.: potentially overrides NC_ERANGE
		 * set by ncmpix_pad_putn_Iint
		 */
		if(lstatus != NC_NOERR)
		{
			ncmpii_free_NC_attr(attrp);
			return lstatus;
		}
	}

	return status;
}

int
ncmpi_get_att_int(int ncid, int varid, const char *name, int *tp)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(attrp->nelems == 0)
		return NC_NOERR;

	if(attrp->type == NC_CHAR)
		return NC_ECHAR;

	{
	const void *xp = attrp->xvalue;
	return ncmpix_pad_getn_Iint(&xp, attrp->nelems, tp, attrp->type);
	}
}


int
ncmpi_put_att_long(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const long *value)
{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;
	NC_attr *old = NULL;
	NC_attr *attrp;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	if(NC_readonly(ncp))
		return NC_EPERM;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	status = ncmpii_cktype(type);
	if(status != NC_NOERR)
		return status;

	if(type == NC_CHAR)
		return NC_ECHAR;

	if((nelems < 0) || (nelems > X_INT_MAX) )/* backward compat */
		return NC_EINVAL; /* Invalid nelems */

	if(nelems != 0 && value == NULL)
		return NC_EINVAL; /* Null arg */

	attrpp = ncmpii_NC_findattr(ncap, name);
	if(attrpp != NULL) /* name in use */
	{
		if(!NC_indef(ncp) )
		{
			const size_t xsz = ncmpix_len_NC_attrV(type, nelems);
			attrp = *attrpp; /* convenience */
	
			if(xsz > attrp->xsz)
				return NC_ENOTINDEFINE;
			/* else, we can reuse existing without redef */
			
			attrp->xsz = xsz;
			attrp->type = type;
			attrp->nelems = nelems;

			if(nelems != 0)
			{
				void *xp = attrp->xvalue;
				status = ncmpix_pad_putn_Ilong(&xp, nelems,
					value, type);
			}
			
			set_NC_hdirty(ncp);

			if(NC_doHsync(ncp))
			{
				const int lstatus = ncmpii_NC_sync(ncp, 1);
				/*
				 * N.B.: potentially overrides NC_ERANGE
				 * set by ncmpix_pad_putn_Ilong
				 */
				if(lstatus != NC_NOERR)
					return lstatus;
			}

			return status;
		}
		/* else, redefine using existing array slot */
		old = *attrpp;
	} 
	else
	{
		if(!NC_indef(ncp))
			return NC_ENOTINDEFINE;

		if(ncap->nelems >= NC_MAX_ATTRS)
			return NC_EMAXATTS;
	}

	status = ncmpii_NC_check_name(name);
	if(status != NC_NOERR)
		return status;

	attrp = ncmpii_new_NC_attr(name, type, nelems);
	if(attrp == NULL)
		return NC_ENOMEM;

	if(nelems != 0)
	{
		void *xp = attrp->xvalue;
		status = ncmpix_pad_putn_Ilong(&xp, nelems,
			value, type);
	}

	if(attrpp != NULL)
	{
		assert(old != NULL);
		*attrpp = attrp;
		ncmpii_free_NC_attr(old);
	}
	else
	{
		const int lstatus = incr_NC_attrarray(ncap, attrp);
		/*
		 * N.B.: potentially overrides NC_ERANGE
		 * set by ncmpix_pad_putn_Ilong
		 */
		if(lstatus != NC_NOERR)
		{
			ncmpii_free_NC_attr(attrp);
			return lstatus;
		}
	}

	return status;
}

int
ncmpi_get_att_long(int ncid, int varid, const char *name, long *tp)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(attrp->nelems == 0)
		return NC_NOERR;

	if(attrp->type == NC_CHAR)
		return NC_ECHAR;

	{
	const void *xp = attrp->xvalue;
	return ncmpix_pad_getn_Ilong(&xp, attrp->nelems, tp, attrp->type);
	}
}


int
ncmpi_put_att_float(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const float *value)
{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;
	NC_attr *old = NULL;
	NC_attr *attrp;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	if(NC_readonly(ncp))
		return NC_EPERM;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	status = ncmpii_cktype(type);
	if(status != NC_NOERR)
		return status;

	if(type == NC_CHAR)
		return NC_ECHAR;

	if((nelems < 0 ) || (nelems > X_INT_MAX)) /* backward compat */
		return NC_EINVAL; /* Invalid nelems */

	if(nelems != 0 && value == NULL)
		return NC_EINVAL; /* Null arg */

	attrpp = ncmpii_NC_findattr(ncap, name);
	if(attrpp != NULL) /* name in use */
	{
		if(!NC_indef(ncp) )
		{
			const size_t xsz = ncmpix_len_NC_attrV(type, nelems);
			attrp = *attrpp; /* convenience */
	
			if(xsz > attrp->xsz)
				return NC_ENOTINDEFINE;
			/* else, we can reuse existing without redef */
			
			attrp->xsz = xsz;
			attrp->type = type;
			attrp->nelems = nelems;

			if(nelems != 0)
			{
				void *xp = attrp->xvalue;
				status = ncmpix_pad_putn_Ifloat(&xp, nelems,
					value, type);
			}
			
			set_NC_hdirty(ncp);

			if(NC_doHsync(ncp))
			{
				const int lstatus = ncmpii_NC_sync(ncp, 1);
				/*
				 * N.B.: potentially overrides NC_ERANGE
				 * set by ncmpix_pad_putn_Ifloat
				 */
				if(lstatus != NC_NOERR)
					return lstatus;
			}

			return status;
		}
		/* else, redefine using existing array slot */
		old = *attrpp;
	} 
	else
	{
		if(!NC_indef(ncp))
			return NC_ENOTINDEFINE;

		if(ncap->nelems >= NC_MAX_ATTRS)
			return NC_EMAXATTS;
	}

	status = ncmpii_NC_check_name(name);
	if(status != NC_NOERR)
		return status;

	attrp = ncmpii_new_NC_attr(name, type, nelems);
	if(attrp == NULL)
		return NC_ENOMEM;

	if(nelems != 0)
	{
		void *xp = attrp->xvalue;
		status = ncmpix_pad_putn_Ifloat(&xp, nelems,
			value, type);
	}

	if(attrpp != NULL)
	{
		assert(old != NULL);
		*attrpp = attrp;
		ncmpii_free_NC_attr(old);
	}
	else
	{
		const int lstatus = incr_NC_attrarray(ncap, attrp);
		/*
		 * N.B.: potentially overrides NC_ERANGE
		 * set by ncmpix_pad_putn_Ifloat
		 */
		if(lstatus != NC_NOERR)
		{
			ncmpii_free_NC_attr(attrp);
			return lstatus;
		}
	}

	return status;
}

int
ncmpi_get_att_float(int ncid, int varid, const char *name, float *tp)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(attrp->nelems == 0)
		return NC_NOERR;

	if(attrp->type == NC_CHAR)
		return NC_ECHAR;

	{
	const void *xp = attrp->xvalue;
	return ncmpix_pad_getn_Ifloat(&xp, attrp->nelems, tp, attrp->type);
	}
}


int
ncmpi_put_att_double(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const double *value)
{
	int status;
	NC *ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;
	NC_attr *old = NULL;
	NC_attr *attrp;

	status = ncmpii_NC_check_id(ncid, &ncp);
	if(status != NC_NOERR)
		return status;

	if(NC_readonly(ncp))
		return NC_EPERM;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	status = ncmpii_cktype(type);
	if(status != NC_NOERR)
		return status;

	if(type == NC_CHAR)
		return NC_ECHAR;

		/* cast needed for braindead systems with signed size_t */
#if 0
	/*commented out because new file format could result in a lot of elements*/
	 if((unsigned long long) nelems > X_INT_MAX) /* backward compat */
        		return NC_EINVAL; /* Invalid nelems */
#endif

	if(nelems != 0 && value == NULL)
		return NC_EINVAL; /* Null arg */

	attrpp = ncmpii_NC_findattr(ncap, name);
	if(attrpp != NULL) /* name in use */
	{
		if(!NC_indef(ncp) )
		{
			const size_t xsz = ncmpix_len_NC_attrV(type, nelems);
			attrp = *attrpp; /* convenience */
	
			if(xsz > attrp->xsz)
				return NC_ENOTINDEFINE;
			/* else, we can reuse existing without redef */
			
			attrp->xsz = xsz;
			attrp->type = type;
			attrp->nelems = nelems;

			if(nelems != 0)
			{
				void *xp = attrp->xvalue;
				status = ncmpix_pad_putn_Idouble(&xp, nelems,
					value, type);
			}
			
			set_NC_hdirty(ncp);

			if(NC_doHsync(ncp))
			{
				const int lstatus = ncmpii_NC_sync(ncp, 1);
				/*
				 * N.B.: potentially overrides NC_ERANGE
				 * set by ncmpix_pad_putn_Idouble
				 */
				if(lstatus != NC_NOERR)
					return lstatus;
			}

			return status;
		}
		/* else, redefine using existing array slot */
		old = *attrpp;
	} 
	else
	{
		if(!NC_indef(ncp))
			return NC_ENOTINDEFINE;

		if(ncap->nelems >= NC_MAX_ATTRS)
			return NC_EMAXATTS;
	}

	status = ncmpii_NC_check_name(name);
	if(status != NC_NOERR)
		return status;

	attrp = ncmpii_new_NC_attr(name, type, nelems);
	if(attrp == NULL)
		return NC_ENOMEM;

	if(nelems != 0)
	{
		void *xp = attrp->xvalue;
		status = ncmpix_pad_putn_Idouble(&xp, nelems,
			value, type);
	}

	if(attrpp != NULL)
	{
		assert(old != NULL);
		*attrpp = attrp;
		ncmpii_free_NC_attr(old);
	}
	else
	{
		const int lstatus = incr_NC_attrarray(ncap, attrp);
		/*
		 * N.B.: potentially overrides NC_ERANGE
		 * set by ncmpix_pad_putn_Idouble
		 */
		if(lstatus != NC_NOERR)
		{
			ncmpii_free_NC_attr(attrp);
			return lstatus;
		}
	}
	return status;
}

int
ncmpi_get_att_double(int ncid, int varid, const char *name, double *tp)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(attrp->nelems == 0)
		return NC_NOERR;

	if(attrp->type == NC_CHAR)
		return NC_ECHAR;

	{
	const void *xp = attrp->xvalue;
	return ncmpix_pad_getn_Idouble(&xp, attrp->nelems, tp, attrp->type);
	}
}



/* deprecated, used to support the 2.x interface */
int
ncmpii_put_att(
	int ncid,
	int varid,
	const char *name,
	nc_type type,
	MPI_Offset nelems,
	const void *value)
{
	switch (type) {
	case NC_BYTE:
		return ncmpi_put_att_schar(ncid, varid, name, type, nelems,
			(schar *)value);
	case NC_CHAR:
		return ncmpi_put_att_text(ncid, varid, name, nelems,
			(char *)value);
	case NC_SHORT:
		return ncmpi_put_att_short(ncid, varid, name, type, nelems,
			(short *)value);
	case NC_INT:
#if (SIZEOF_INT >= X_SIZEOF_INT)
		return ncmpi_put_att_int(ncid, varid, name, type, nelems,
			(int *)value);
#elif SIZEOF_LONG == X_SIZEOF_INT
		return ncmpi_put_att_long(ncid, varid, name, type, nelems,
			(long *)value);
#endif
	case NC_FLOAT:
		return ncmpi_put_att_float(ncid, varid, name, type, nelems,
			(float *)value);
	case NC_DOUBLE:
		return ncmpi_put_att_double(ncid, varid, name, type, nelems,
			(double *)value);
	default:
		return NC_EBADTYPE;
	}
}


/* deprecated, used to support the 2.x interface */
int
ncmpii_get_att(int ncid, int varid, const char *name, void *value)
{
	int status;
	NC_attr *attrp;

	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	switch (attrp->type) {
	case NC_BYTE:
		return ncmpi_get_att_schar(ncid, varid, name,
			(schar *)value);
	case NC_CHAR:
		return ncmpi_get_att_text(ncid, varid, name,
			(char *)value);
	case NC_SHORT:
		return ncmpi_get_att_short(ncid, varid, name,
			(short *)value);
	case NC_INT:
#if (SIZEOF_INT >= X_SIZEOF_INT)
		return ncmpi_get_att_int(ncid, varid, name,
			(int *)value);
#elif SIZEOF_LONG == X_SIZEOF_INT
		return ncmpi_get_att_long(ncid, varid, name,
			(long *)value);
#endif
	case NC_FLOAT:
		return ncmpi_get_att_float(ncid, varid, name,
			(float *)value);
	case NC_DOUBLE:
		return ncmpi_get_att_double(ncid, varid, name,
			(double *)value);
	default:
		return NC_EBADTYPE;
	}
}
