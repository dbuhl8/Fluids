/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: var.c 829 2010-05-26 20:17:57Z wkliao $ */

#include "nc.h"
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h>
#include <assert.h>
#include "ncx.h"
#include "rnd.h"
#include "macro.h"

/* Prototypes for functions used only in this file */
static MPI_Offset ncx_szof(nc_type type);

/*
 * Free var
 * Formerly
NC_free_var(var)
 */
void
ncmpii_free_NC_var(NC_var *varp)
{
	if(varp == NULL)
		return;
	ncmpii_free_NC_attrarrayV(&varp->attrs);
	ncmpii_free_NC_string(varp->name);
	NCI_Free(varp);
}


/* 
 * Common code for ncmpii_new_NC_var() 
 * and ncx_get_NC_var()
 */
NC_var *
ncmpii_new_x_NC_var(
	NC_string *strp,
	size_t ndims)
{
	NC_var *varp;
	const MPI_Offset o1 = M_RNDUP(ndims * sizeof(MPI_Offset));
	const MPI_Offset o2 = M_RNDUP(ndims * sizeof(MPI_Offset));
	const MPI_Offset sz =  M_RNDUP(sizeof(NC_var)) +
		 o1 + o2 + ndims * sizeof(MPI_Offset);

	varp = (NC_var *) NCI_Malloc(sz);
	if(varp == NULL )
		return NULL;
	(void) memset(varp, 0, sz);

	varp->name = strp;
	varp->ndims = ndims;

	if(ndims != 0)
	{
		/*
		 * NOTE: lint may complain about the next 3 lines:
		 * "pointer cast may result in improper alignment".
		 * We use the M_RNDUP() macro to get the proper alignment.
		 */
		varp->dimids = (int *)((char *)varp + M_RNDUP(sizeof(NC_var)));
		varp->shape = (MPI_Offset *)((char *)varp->dimids + o1);
		varp->dsizes = (MPI_Offset *)((char *)varp->shape + o2);
	}


	varp->xsz = 0;
	varp->len = 0;
	varp->begin = 0;

	return varp;
}


/*
 * Formerly
NC_new_var()
 */
static NC_var *
ncmpii_new_NC_var(const char *name, nc_type type,
	size_t ndims, const int *dimids)
{
	NC_string *strp;
	NC_var *varp;

	strp = ncmpii_new_NC_string(strlen(name), name);
	if(strp == NULL)
		return NULL;

	varp = ncmpii_new_x_NC_var(strp, ndims);
	if(varp == NULL )
	{
		ncmpii_free_NC_string(strp);
		return NULL;
	}
	
	varp->type = type;

	if( ndims != 0 && dimids != NULL)
		(void) memcpy(varp->dimids, dimids, ndims * sizeof(int));

	return(varp);
}

static NC_var *
dup_NC_var(const NC_var *rvarp)
{
	NC_var *varp = ncmpii_new_NC_var(rvarp->name->cp, rvarp->type,
		 rvarp->ndims, rvarp->dimids);
	if(varp == NULL)
		return NULL;

	
	if(ncmpii_dup_NC_attrarrayV(&varp->attrs, &rvarp->attrs) != NC_NOERR)
	{
		ncmpii_free_NC_var(varp);
		return NULL;
	}

	(void) memcpy(varp->shape, rvarp->shape,
			 rvarp->ndims * sizeof(MPI_Offset));
	(void) memcpy(varp->dsizes, rvarp->dsizes,
			 rvarp->ndims * sizeof(MPI_Offset));
	varp->xsz = rvarp->xsz;
	varp->len = rvarp->len;
	varp->begin = rvarp->begin;

	return varp;
}


/* vararray */


/*
 * Free the stuff "in" (referred to by) an NC_vararray.
 * Leaves the array itself allocated.
 */
void
ncmpii_free_NC_vararrayV0(NC_vararray *ncap)
{
	assert(ncap != NULL);

	if(ncap->nelems == 0)
		return;

	assert(ncap->value != NULL);

	{
		NC_var **vpp = ncap->value;
		NC_var *const *const end = &vpp[ncap->nelems];
		for( /*NADA*/; vpp < end; vpp++)
		{
			ncmpii_free_NC_var(*vpp);
			*vpp = NULL;
		}
	}
	ncap->nelems = 0;
}


/*
 * Free NC_vararray values.
 * formerly
NC_free_array()
 */
void
ncmpii_free_NC_vararrayV(NC_vararray *ncap)
{
	assert(ncap != NULL);
	
	if(ncap->nalloc == 0)
		return;

	assert(ncap->value != NULL);

	ncmpii_free_NC_vararrayV0(ncap);

	NCI_Free(ncap->value);
	ncap->value = NULL;
	ncap->nalloc = 0;
}


int
ncmpii_dup_NC_vararrayV(NC_vararray *ncap, const NC_vararray *ref)
{
	int status = NC_NOERR;

	assert(ref != NULL);
	assert(ncap != NULL);

	if(ref->nelems != 0)
	{
		const MPI_Offset sz = ref->nelems * sizeof(NC_var *);
		ncap->value = (NC_var **) NCI_Malloc(sz);
		if(ncap->value == NULL)
			return NC_ENOMEM;
		(void) memset(ncap->value, 0, sz);
		ncap->nalloc = ref->nelems;
	}

	ncap->nelems = 0;
	{
		NC_var **vpp = ncap->value;
		const NC_var **drpp = (const NC_var **)ref->value;
		NC_var *const *const end = &vpp[ref->nelems];
		for( /*NADA*/; vpp < end; drpp++, vpp++, ncap->nelems++)
		{
			*vpp = dup_NC_var(*drpp);
			if(*vpp == NULL)
			{
				status = NC_ENOMEM;
				break;
			}
		}
	}

	if(status != NC_NOERR)
	{
		ncmpii_free_NC_vararrayV(ncap);
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
incr_NC_vararray(NC_vararray *ncap, NC_var *newelemp)
{
	NC_var **vp;

	assert(ncap != NULL);

	if(ncap->nalloc == 0)
	{
		assert(ncap->nelems == 0);
		vp = (NC_var **) NCI_Malloc(NC_ARRAY_GROWBY * sizeof(NC_var *));
		if(vp == NULL)
			return NC_ENOMEM;
		ncap->value = vp;
		ncap->nalloc = NC_ARRAY_GROWBY;
	}
	else if(ncap->nelems +1 > ncap->nalloc)
	{
		vp = (NC_var **) NCI_Realloc(ncap->value,
			(ncap->nalloc + NC_ARRAY_GROWBY) * sizeof(NC_var *));
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


static NC_var *
elem_NC_vararray(const NC_vararray *ncap, MPI_Offset elem)
{
	assert(ncap != NULL);
		/* cast needed for braindead systems with signed MPI_Offset */
	if((elem < 0) ||  ncap->nelems == 0 || elem >= ncap->nelems)
		return NULL;

	assert(ncap->value != NULL);

	return ncap->value[elem];
}


/* End vararray per se */


/*
 * Step thru NC_VARIABLE array, seeking match on name.
 * Return varid or -1 on not found.
 * *varpp is set to the appropriate NC_var.
 * Formerly (sort of)
NC_hvarid
 */
int
ncmpii_NC_findvar(const NC_vararray *ncap, const char *name, NC_var **varpp)
{
	NC_var **loc;
	MPI_Offset slen;
	int varid;

	assert(ncap != NULL);

	if(ncap->nelems == 0)
		return -1;

	loc = (NC_var **) ncap->value;

	slen = strlen(name);

	for(varid = 0; (MPI_Offset) varid < ncap->nelems; varid++, loc++)
	{
		if(strlen((*loc)->name->cp) == slen &&
			strncmp((*loc)->name->cp, name, slen) == 0)
		{
			if(varpp != NULL)
				*varpp = *loc;
			return(varid); /* Normal return */
		}
	}
	return(-1); /* not found */
}

/* 
 * For a netcdf type
 *  return the size of one element in the external representation.
 * Note that arrays get rounded up to X_ALIGN boundaries.
 * Formerly
NC_xtypelen
 * See also ncx_len()
 */
static MPI_Offset
ncx_szof(nc_type type)
{
	switch(type){
	case NC_BYTE:
	case NC_CHAR:
		return(1);
	case NC_SHORT :
		return(2);
	case NC_INT:
		return X_SIZEOF_INT;
	case NC_FLOAT:
		return X_SIZEOF_FLOAT;
	case NC_DOUBLE : 
		return X_SIZEOF_DOUBLE;
	default:
		assert("ncx_szof invalid type" == 0);
	}
	/* default */
	return 0;
}


/*
 * 'compile' the shape and len of a variable
 *  Formerly
 *   ncmpii_NC_var_shape(var, dims)
 */
#include <stdio.h>
int
ncmpii_NC_var_shape64(NC_var *varp, const NC_dimarray *dims)
{
	MPI_Offset *shp, *dsp, *op;
	int *ip;
	const NC_dim *dimp;
	MPI_Offset product = 1;
/*
	int *shp, *dsp, *op;
	int *ip;
	const NC_dim *dimp;
	int product = 1;
*/	
	varp->xsz = ncx_szof(varp->type);

	if(varp->ndims == 0)
	{
		goto out;
	}

	/*
	 * use the user supplied dimension indices
	 * to determine the shape
	 */
/* 	for (i=0; i<varp->ndims;i++){
  		   
                   printf("ncmpii_NC_var_shape64: varp->dimids[%d]:%d, varp->shape[%d]:%d\n",i, varp->dimids[i],i, varp->shape[i]);

        };
*/
        for(ip = varp->dimids, op = varp->shape
		; ip < &varp->dimids[varp->ndims]; ip++, op++)
	{
		if(*ip < 0 || (*ip) >= ((dims != NULL) ? dims->nelems : 1) ){
#ifdef DEBUG
		        printf("return NC_EBADDIM, *ip:%d\n", *ip);
#endif
			return NC_EBADDIM;
		}
		dimp = ncmpii_elem_NC_dimarray(dims, *ip);
		*op = dimp->size;
		if(*op == NC_UNLIMITED && ip != varp->dimids)
			return NC_EUNLIMPOS;
	}

	/* 
	 * Compute the dsizes
	 */
				/* ndims is > 0 here */
	for(shp = varp->shape + varp->ndims -1,
				dsp = varp->dsizes + varp->ndims -1;
 			shp >= varp->shape;
			shp--, dsp--)
	{
		if(!(shp == varp->shape && IS_RECVAR(varp)))
		{
			if( *shp <= X_UINT_MAX / product)
			{
				product *= *shp;
			} else
			{
				product = X_UINT_MAX;
			}
		}
		*dsp = product;
	}


out :

	/* with new CDF-5 variable sizes, we might be able to avoid the bug
	 * from serial netcdf.  not sure. need to double check */
	if (varp->xsz <= X_UINT_MAX - 1/ product) 
		/* if int. mult won't overflow ...*/
	{
		varp->len = product * varp->xsz;
		switch(varp->type) {
			case NC_BYTE :
			case NC_CHAR : 
			case NC_SHORT :
				if( varp->len%4 != 0 )
				{
					varp->len += 4 - varp->len%4; /* round up */
					/*		*dsp += 4 - *dsp%4; */
				}
				break;
			default:
				/* already aligned */
				break;
		}
	} else
	{ /* ok for last var to be "too big", indicated by this special len */
		varp->len = X_UINT_MAX;
	}
#if 0
	arrayp("\tshape", varp->ndims, varp->shape);
	arrayp("\tdsizes", varp->ndims, varp->dsizes);
#endif
	return NC_NOERR;
}

/*
 * Check whether variable size is less than or equal to vlen_max,
 * without overflowing in arithmetic calculations.  If OK, return 1,
 * else, return 0.  For CDF1 format or for CDF2 format on non-LFS
 * platforms, vlen_max should be 2^31 - 4, but for CDF2 format on
 * systems with LFS it should be 2^32 - 4.
 */
int
ncmpii_NC_check_vlen(NC_var *varp, MPI_Offset vlen_max) {
    MPI_Offset prod=varp->xsz;     /* product of xsz and dimensions so far */

    int ii;

    for(ii = IS_RECVAR(varp) ? 1 : 0; ii < varp->ndims; ii++) {
       if (varp->shape[ii] > vlen_max / prod) {
           return 0;           /* size in bytes won't fit in a 32-bit int */
       }
       prod *= varp->shape[ii];
    }
    return 1;                  /* OK */
}

/*
 * Given valid ncp and varid, return var
 *  else NULL on error
 * Formerly
NC_hlookupvar()
 */
NC_var *
ncmpii_NC_lookupvar(NC *ncp,  int varid)
{
	NC_var *varp;

	if(varid == NC_GLOBAL)
	{
		/* Global is error in this context */
		return(NULL);
	}

	varp = elem_NC_vararray(&ncp->vars, varid);
	if(varp == NULL)
	{
		return NULL;
	}

	assert(varp != NULL);

	return(varp);
}


/* Public */

int
ncmpi_def_var( int ncid, const char *name, nc_type type,
	 int ndims, const int *dimids, int *varidp)
{
	int status;
	NC *ncp;
	int varid;
	NC_var *varp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;


	if(!NC_indef(ncp))
	{
		return NC_ENOTINDEFINE;
	}

	status = ncmpii_NC_check_name(name);
	if(status != NC_NOERR)
		return status;

	status = ncmpii_cktype(type);
	if(status != NC_NOERR)
		return status;

	if((ndims < 0) || ndims > X_INT_MAX) /* Backward compat */
	{
		return NC_EINVAL;
	} 

	if(ncp->vars.nelems >= NC_MAX_VARS)
	{
		return NC_EMAXVARS;
	}

	varid = ncmpii_NC_findvar(&ncp->vars, name, &varp);
	if(varid != -1)
	{
		return NC_ENAMEINUSE;
	}
	

	
	varp = ncmpii_new_NC_var(name, type, ndims, dimids);
	if(varp == NULL)
		return NC_ENOMEM;

	

	status = ncmpii_NC_var_shape64(varp, &ncp->dims);
	if(status != NC_NOERR)
	{
		ncmpii_free_NC_var(varp);
		return status;
	}

	status = incr_NC_vararray(&ncp->vars, varp);
	if(status != NC_NOERR)
	{
		ncmpii_free_NC_var(varp);
		return status;
	}

	if(varidp != NULL)
		*varidp = (int)ncp->vars.nelems -1; /* varid */
	return NC_NOERR;
}


int
ncmpi_inq_varid(int ncid, const char *name, int *varid_ptr)
{
	int status;
	NC *ncp;
	NC_var *varp;
	int varid;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	varid = ncmpii_NC_findvar(&ncp->vars, name, &varp);
	if(varid == -1)
	{
		return NC_ENOTVAR;
	}

	*varid_ptr = varid;
	return NC_NOERR;
}


int
ncmpi_inq_var(int ncid,
	int varid,
	char *name,
	nc_type *typep,
	int *ndimsp,
	int *dimids,
	int *nattsp)
{
	int status;
	NC *ncp;
	NC_var *varp;
	MPI_Offset ii;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	varp = elem_NC_vararray(&ncp->vars, varid);
	if(varp == NULL)
		return NC_ENOTVAR;

	if(name != NULL)
	{
		(void) strncpy(name, varp->name->cp, varp->name->nchars);
		name[varp->name->nchars] = 0;
	}

	if(typep != 0)
		*typep = varp->type;
	if(ndimsp != 0)
	{
		*ndimsp = varp->ndims;
	}
	if(dimids != 0)
	{
		for(ii = 0; ii < varp->ndims; ii++)
		{
			dimids[ii] = varp->dimids[ii];
		}
	}
	if(nattsp != 0)
	{
		*nattsp = (int) varp->attrs.nelems;
	}

	return NC_NOERR;
}


int 
ncmpi_inq_varname(int ncid,  int varid, char *name)
{
	int status;
	NC *ncp;
	NC_var *varp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	varp = elem_NC_vararray(&ncp->vars, varid);
	if(varp == NULL)
		return NC_ENOTVAR;

	if(name != NULL)
	{
		(void) strncpy(name, varp->name->cp, varp->name->nchars);
		name[varp->name->nchars] = 0;
	}

	return NC_NOERR;
}

int 
ncmpi_inq_vartype(int ncid,  int varid, nc_type *typep)
{
	int status;
	NC *ncp;
	NC_var *varp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	varp = elem_NC_vararray(&ncp->vars, varid);
	if(varp == NULL)
		return NC_ENOTVAR;

	if(typep != 0)
		*typep = varp->type;

	return NC_NOERR;
}

int 
ncmpi_inq_varndims(int ncid,  int varid, int *ndimsp)
{
	int status;
	NC *ncp;
	NC_var *varp;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	varp = elem_NC_vararray(&ncp->vars, varid);
	if(varp == NULL)
		return NC_ENOTVAR; /* TODO: is this the right error code? */

	if(ndimsp != 0)
	{
		*ndimsp = (int) varp->ndims;
	}

	return NC_NOERR;
}


int 
ncmpi_inq_vardimid(int ncid,  int varid, int *dimids)
{
	int status;
	NC *ncp;
	NC_var *varp;
	MPI_Offset ii;

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	varp = elem_NC_vararray(&ncp->vars, varid);
	if(varp == NULL)
		return NC_ENOTVAR; /* TODO: is this the right error code? */

	if(dimids != 0)
	{
		for(ii = 0; ii < varp->ndims; ii++)
		{
			dimids[ii] = varp->dimids[ii];
		}
	}

	return NC_NOERR;
}


int 
ncmpi_inq_varnatts(int ncid,  int varid, int *nattsp)
{
	int status;
	NC *ncp;
	NC_var *varp;

	if(varid == NC_GLOBAL)
		return ncmpi_inq_natts(ncid, nattsp);

	status = ncmpii_NC_check_id(ncid, &ncp); 
	if(status != NC_NOERR)
		return status;

	varp = elem_NC_vararray(&ncp->vars, varid);
	if(varp == NULL)
		return NC_ENOTVAR; /* TODO: is this the right error code? */

	if(nattsp != 0)
	{
		*nattsp = (int) varp->attrs.nelems;
	}

	return NC_NOERR;
}

int
ncmpi_rename_var(int ncid,  int varid, const char *newname)
{
    int status, other;
    NC *ncp;
    NC_var *varp;
    NC_string *old, *newStr;

    status = ncmpii_NC_check_id(ncid, &ncp); 
    if (status != NC_NOERR)
        return status;

    if (NC_readonly(ncp))
        return NC_EPERM;

    status = ncmpii_NC_check_name(newname);
    if (status != NC_NOERR)
        return status;

    /* check for name in use */
    other = ncmpii_NC_findvar(&ncp->vars, newname, &varp);
    if (other != -1)
        return NC_ENAMEINUSE;
        
    varp = ncmpii_NC_lookupvar(ncp, varid);
    if (varp == NULL)
        /* invalid varid */
        return NC_ENOTVAR; /* TODO: is this the right error code? */

    old = varp->name;
    if (NC_indef(ncp)) {
       newStr = ncmpii_new_NC_string(strlen(newname),newname);
       if (newStr == NULL)
           return(-1);
       varp->name = newStr;
       ncmpii_free_NC_string(old);
       return NC_NOERR;
    }
    /* else, not in define mode */

    status = ncmpii_set_NC_string(varp->name, newname);
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
