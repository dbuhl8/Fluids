/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: nc.h 839 2010-06-03 17:49:17Z robl $ */
#ifndef _NC_H_
#define _NC_H_

/*
 *	netcdf library 'private' data structures, objects and interfaces
 */

#include "ncconfig.h"

#include <stddef.h>	/* MPI_Offset */
#include <sys/types.h>	/* MPI_Offset */
#include "pnetcdf.h"
#include "ncio.h"	/* ncio */
#include "fbits.h"


/* XXX: this seems really low.  do we end up spending a ton of time mallocing?
 * could we reduce that by increasing this to something 21st century? */
#ifndef NC_ARRAY_GROWBY
#define NC_ARRAY_GROWBY 4
#endif

/* ncmpi_create/ncmpi_open set up header to be 'chunksize' big and to grow by 'chunksize' as new items added */
#define NC_DEFAULT_CHUNKSIZE 4096

/*
 * The extern size of an empty
 * netcdf version 1 file.
 * The initial value of ncp->xsz.
 */
#define MIN_NC_XSZ 32

typedef struct NC NC; /* forward reference */

/*
 *  The internal data types
 */
typedef enum {
	NC_UNSPECIFIED = 0,
/* future	NC_BITFIELD = 7, */
/*	NC_STRING =	8,	*/
	NC_DIMENSION =	10,
	NC_VARIABLE =	11,
	NC_ATTRIBUTE =	12
} NCtype;


/*
 * Counted string for names and such
 */
typedef struct {
	/* all xdr'd */
	MPI_Offset nchars;
	char *cp;
} NC_string;

extern NC *
ncmpii_new_NC(const MPI_Offset *chunkp);

extern NC *
ncmpii_dup_NC(const NC *ref);

/* Begin defined in string.c */
extern void
ncmpii_free_NC_string(NC_string *ncstrp);

extern int
ncmpii_NC_check_name(const char *name);

extern NC_string *
ncmpii_new_NC_string(MPI_Offset slen, const char *str);

extern int
ncmpii_set_NC_string(NC_string *ncstrp, const char *str);

/* End defined in string.c */

/*
 * NC dimension stucture
 */
typedef struct {
	/* all xdr'd */
	NC_string *name;
	MPI_Offset size;
} NC_dim;

typedef struct NC_dimarray {
	size_t nalloc;		/* number allocated >= nelems */
	/* below gets xdr'd */
	/* NCtype type = NC_DIMENSION */
	MPI_Offset nelems;		/* number of defined variables */
	NC_dim **value;
} NC_dimarray;

/* Begin defined in dim.c */

extern void
ncmpii_free_NC_dim(NC_dim *dimp);

extern NC_dim *
ncmpii_new_x_NC_dim(NC_string *name);

extern int
ncmpii_find_NC_Udim(const NC_dimarray *ncap, NC_dim **dimpp);

/* dimarray */

extern void
ncmpii_free_NC_dimarrayV0(NC_dimarray *ncap);

extern void
ncmpii_free_NC_dimarrayV(NC_dimarray *ncap);

extern int
ncmpii_dup_NC_dimarrayV(NC_dimarray *ncap, const NC_dimarray *ref);

extern NC_dim *
ncmpii_elem_NC_dimarray(const NC_dimarray *ncap, size_t elem);

extern int
ncmpi_def_dim(int ncid, const char *name, MPI_Offset size, int *dimidp);

extern int
ncmpi_rename_dim( int ncid, int dimid, const char *newname);

extern int
ncmpi_inq_dimid(int ncid, const char *name, int *dimid_ptr);

extern int
ncmpi_inq_dim(int ncid, int dimid, char *name, MPI_Offset *sizep);

extern int 
ncmpi_inq_dimname(int ncid, int dimid, char *name);

extern int 
ncmpi_inq_dimlen(int ncid, int dimid, MPI_Offset *lenp);
/* End defined in dim.c */

/*
 * NC attribute
 */
typedef struct {
	MPI_Offset xsz;		/* amount of space at xvalue */
	/* below gets xdr'd */
	NC_string *name;
	nc_type type;		/* the discriminant */
	MPI_Offset nelems;	/* number of defined variables */
	void *xvalue;		/* the actual data, in external representation */
} NC_attr;

typedef struct NC_attrarray {
	MPI_Offset nalloc;	/* number allocated >= nelems */
	/* below gets xdr'd */
	/* NCtype type = NC_ATTRIBUTE */
	MPI_Offset nelems;	/* number of defined variables */
	NC_attr **value;
} NC_attrarray;

/* Begin defined in attr.c */

extern void
ncmpii_free_NC_attr(NC_attr *attrp);

extern NC_attr *
ncmpii_new_x_NC_attr(
	NC_string *strp,
	nc_type type,
	MPI_Offset nelems);

extern NC_attr **
ncmpii_NC_findattr(const NC_attrarray *ncap, const char *name);

/* attrarray */

extern void
ncmpii_free_NC_attrarrayV0(NC_attrarray *ncap);

extern void
ncmpii_free_NC_attrarrayV(NC_attrarray *ncap);

extern int
ncmpii_dup_NC_attrarrayV(NC_attrarray *ncap, const NC_attrarray *ref);

extern NC_attr *
ncmpii_elem_NC_attrarray(const NC_attrarray *ncap, MPI_Offset elem);

extern int
ncmpi_put_att_text(int ncid, int varid, const char *name,
	MPI_Offset nelems, const char *value);

extern int
ncmpi_get_att_text(int ncid, int varid, const char *name, char *str);

extern int
ncmpi_put_att_schar(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const signed char *value);

extern int
ncmpi_get_att_schar(int ncid, int varid, const char *name, signed char *tp);

extern int
ncmpi_put_att_uchar(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const unsigned char *value);

extern int
ncmpi_get_att_uchar(int ncid, int varid, const char *name, unsigned char *tp);

extern int
ncmpi_put_att_short(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const short *value);

extern int
ncmpi_get_att_short(int ncid, int varid, const char *name, short *tp);

extern int
ncmpi_put_att_int(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const int *value);

extern int
ncmpi_get_att_int(int ncid, int varid, const char *name, int *tp);

extern int
ncmpi_put_att_long(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const long *value);

extern int
ncmpi_get_att_long(int ncid, int varid, const char *name, long *tp);

extern int
ncmpi_put_att_float(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const float *value);
extern int
ncmpi_get_att_float(int ncid, int varid, const char *name, float *tp);
extern int
ncmpi_put_att_double(int ncid, int varid, const char *name,
	nc_type type, MPI_Offset nelems, const double *value);
extern int
ncmpi_get_att_double(int ncid, int varid, const char *name, double *tp);

extern int 
ncmpi_inq_attid(int ncid, int varid, const char *name, int *attnump);

extern int 
ncmpi_inq_atttype(int ncid, int varid, const char *name, nc_type *datatypep);

extern int 
ncmpi_inq_attlen(int ncid, int varid, const char *name, MPI_Offset *lenp);

extern int
ncmpi_inq_att(int ncid, int varid, const char *name, 
	nc_type *datatypep, MPI_Offset *lenp);

extern int
ncmpi_copy_att(int ncid_in, int varid_in, const char *name, 
		int ncid_out, int ovarid);

extern int
ncmpi_rename_att( int ncid, int varid, const char *name, const char *newname);

extern int
ncmpi_del_att(int ncid, int varid, const char *name);

extern int
ncmpi_inq_attname(int ncid, int varid, int attnum, char *name);
/* End defined in attr.c */
/*
 * NC variable: description and data
 */
typedef struct {
	MPI_Offset xsz;		/* xszof 1 element */
	MPI_Offset *shape; /* compiled info: dim->size of each dim */
	MPI_Offset *dsizes; /* compiled info: the right to left product of shape */
	/* below gets xdr'd */
	NC_string *name;
	/* next two: formerly NC_iarray *assoc */ /* user definition */
	size_t ndims;	/* assoc->count */
	int *dimids;	/* assoc->value */
	NC_attrarray attrs;
	nc_type type;		/* the discriminant */
	MPI_Offset len;		/* the total length originally allocated */
	MPI_Offset begin;
} NC_var;

typedef struct NC_vararray {
	MPI_Offset nalloc;	/* number allocated >= nelems */
	/* below gets xdr'd */
	/* NCtype type = NC_VARIABLE */
	MPI_Offset nelems;	/* number of defined variables */
	NC_var **value;
} NC_vararray;

/* Begin defined in var.c */

extern void
ncmpii_free_NC_var(NC_var *varp);

extern NC_var *
ncmpii_new_x_NC_var(
	NC_string *strp,
	size_t ndims
        );

/* vararray */

extern void
ncmpii_free_NC_vararrayV0(NC_vararray *ncap);

extern void
ncmpii_free_NC_vararrayV(NC_vararray *ncap);

extern int
ncmpii_dup_NC_vararrayV(NC_vararray *ncap, const NC_vararray *ref);

extern int
ncmpii_NC_var_shape64(NC_var *varp, const NC_dimarray *dims);

extern int
ncmpii_NC_findvar(const NC_vararray *ncap, const char *name, NC_var **varpp);

extern int
ncmpii_NC_check_vlen(NC_var *varp, MPI_Offset vlen_max);

extern NC_var *
ncmpii_NC_lookupvar(NC *ncp, int varid);

extern int
ncmpi_def_var( int ncid, const char *name, nc_type type,
              int ndims, const int *dimidsp, int *varidp);

extern int
ncmpi_rename_var(int ncid, int varid, const char *newname);

extern int
ncmpi_inq_var(int ncid, int varid, char *name, nc_type *typep, 
		int *ndimsp, int *dimids, int *nattsp);

extern int
ncmpi_inq_varid(int ncid, const char *name, int *varid_ptr);

extern int 
ncmpi_inq_varname(int ncid, int varid, char *name);

extern int 
ncmpi_inq_vartype(int ncid, int varid, nc_type *typep);

extern int 
ncmpi_inq_varndims(int ncid, int varid, int *ndimsp);

extern int 
ncmpi_inq_vardimid(int ncid, int varid, int *dimids);

extern int 
ncmpi_inq_varnatts(int ncid, int varid, int *nattsp);

extern int
ncmpi_rename_var(int ncid, int varid, const char *newname);
/* End defined in var.c */

#define IS_RECVAR(vp) \
	((vp)->shape != NULL ? (*(vp)->shape == NC_UNLIMITED) : 0 )

/*
 *  *  The PnetCDF non-blocking I/O request type
 *   */
typedef struct NC_req {
    int            id;
    int            rw_flag;
    NC_var        *varp;
    void          *buf;
    void          *xbuf;
    void          *cbuf;
    void          *lbuf;
    int            iscontig_of_ptypes;
    int            is_imap;
    int            ndims;
    MPI_Offset    *start;  /* [ndims] */
    MPI_Offset    *count;  /* [ndims] */
    MPI_Offset    *stride; /* [ndims] */
    MPI_Offset     nelems;
    MPI_Offset     cnelems;
    MPI_Offset     lnelems;
    MPI_Offset     bufcount;
    MPI_Offset     offset_start;  /* starting of aggregate access region */
    MPI_Offset     offset_end;    /*   ending of aggregate access region */
    MPI_Datatype   datatype;
    MPI_Datatype   ptype;
    MPI_Datatype   imaptype;
    int           *status;
    int            num_subreqs;
    struct NC_req *subreqs;  /* [num_subreq] */
    struct NC_req *next;
} NC_req;

struct NC {
	/* links to make list of open netcdf's */
	struct NC *next;
	struct NC *prev;
	/* contains the previous NC during redef. */
	struct NC *old;
	/* flags */
#define NC_INDEP 1	/* in independent data mode, cleared by endindep */
#define NC_CREAT 2	/* in create phase, cleared by ncenddef */
#define NC_INDEF 8	/* in define mode, cleared by ncenddef */
#define NC_NSYNC 0x10	/* synchronise numrecs on change */
#define NC_HSYNC 0x20	/* synchronise whole header on change */
#define NC_NDIRTY 0x40	/* numrecs has changed */
#define NC_HDIRTY 0x80  /* header info has changed */
/*	NC_NOFILL in netcdf.h, historical interface */
	int flags;
	ncio *nciop;
	MPI_Offset chunk;	/* largest extent this layer will request from ncio->get() */
	MPI_Offset xsz;	/* external size of this header, <= var[0].begin */
	MPI_Offset begin_var; /* position of the first (non-record) var */
	MPI_Offset begin_rec; /* position of the first 'record' */
	/* don't constrain maximum size of record unnecessarily */
	MPI_Offset recsize;	/* length of 'record': sum of single record sizes from all record variables */	
	/* below gets xdr'd */
	MPI_Offset numrecs; /* number of 'records' allocated */
	NC_dimarray dims;
	NC_attrarray attrs;
	NC_vararray vars;
        NC_req *head;
        NC_req *tail;
};

#define NC_readonly(ncp) \
	(!fIsSet(ncp->nciop->ioflags, NC_WRITE))

#define NC_IsNew(ncp) \
	fIsSet((ncp)->flags, NC_CREAT)

#define NC_indep(ncp) \
	fIsSet((ncp)->flags, NC_INDEP)

#define NC_indef(ncp) \
	(NC_IsNew(ncp) || fIsSet((ncp)->flags, NC_INDEF)) 

#define set_NC_ndirty(ncp) \
	fSet((ncp)->flags, NC_NDIRTY)

#define NC_ndirty(ncp) \
	fIsSet((ncp)->flags, NC_NDIRTY)

#define set_NC_hdirty(ncp) \
	fSet((ncp)->flags, NC_HDIRTY)

#define NC_hdirty(ncp) \
	fIsSet((ncp)->flags, NC_HDIRTY)

#define NC_dofill(ncp) \
	(!fIsSet((ncp)->flags, NC_NOFILL))

#define NC_doHsync(ncp) \
	fIsSet((ncp)->flags, NC_HSYNC)

#define NC_doNsync(ncp) \
	fIsSet((ncp)->flags, NC_NSYNC)

#define NC_get_numrecs(ncp) \
	((ncp)->numrecs)

#define NC_set_numrecs(ncp, nrecs) \
	{((ncp)->numrecs = (nrecs));}

#define NC_increase_numrecs(ncp, nrecs) \
	{if((nrecs) > (ncp)->numrecs) ((ncp)->numrecs = (nrecs));}
/* Begin defined in nc.c */

extern int
ncmpii_NC_check_id(int ncid, NC **ncpp);

extern int
ncmpii_cktype(nc_type datatype);

extern MPI_Offset
ncmpix_howmany(nc_type type, MPI_Offset xbufsize);

extern int
ncmpii_read_numrecs(NC *ncp);

extern int
ncmpii_write_numrecs(NC *ncp);

extern int
ncmpii_NC_sync(NC *ncp, int doFsync);

extern void
ncmpii_free_NC(NC *ncp);

extern void
ncmpii_add_to_NCList(NC *ncp);

extern void
ncmpii_del_from_NCList(NC *ncp);

extern int
ncmpii_read_NC(NC *ncp);

extern int 
ncmpii_NC_enddef(NC *ncp);

extern int 
ncmpii_NC_close(NC *ncp);

extern int
ncmpi_inq(int ncid, int *ndimsp, int *nvarsp, int *nattsp, int *xtendimp);

extern int 
ncmpi_inq_ndims(int ncid, int *ndimsp);

extern int 
ncmpi_inq_nvars(int ncid, int *nvarsp);

extern int 
ncmpi_inq_natts(int ncid, int *nattsp);

extern int 
ncmpi_inq_unlimdim(int ncid, int *xtendimp);

extern int
ncmpi_get_default_format(void);

/* End defined in nc.c */
/* Begin defined in v1hpg.c */

extern size_t
ncx_len_NC(const NC *ncp, MPI_Offset sizeof_off_t);

extern int
ncx_put_NC(const NC *ncp, void **xpp, MPI_Offset offset, MPI_Offset extent);

extern int
nc_get_NC( NC *ncp);

/* End defined in v1hpg.c */

#if 0
/* Begin defined in putget.c */

extern int
ncmpii_fill_NC_var(NC *ncp, const NC_var *varp, MPI_Offset recno);

extern int
ncmpii_inq_rec(int ncid, MPI_Offset *nrecvars, MPI_Offset *recvarids, MPI_Offset *recsizes);

extern int
ncmpii_get_rec(int ncid, MPI_Offset recnum, void **datap);

extern int
ncmpii_put_rec(int ncid, MPI_Offset recnum, void *const *datap);
#endif

/* End defined in putget.c */

/* Begin defined in header.c */
typedef struct bufferinfo {
  ncio *nciop;		
  MPI_Offset offset;	/* current read/write offset in the file */
  int version;		/* either 1 for normal netcdf or 
			   2 for 8-byte offset version, 
			   5 for NC_FORMAT_64BIT_DATA version */
  void *base;     	/* beginning of read/write buffer */
  void *pos;      	/* current position in buffer */
  MPI_Offset size;		/* size of the buffer */
  MPI_Offset index;		/* index of current position in buffer */
} bufferinfo;  

extern MPI_Offset 
ncmpix_len_nctype(nc_type type);

#if 0
extern int
hdr_put_NC_attrarray(bufferinfo *pbp, const NC_attrarray *ncap);
#endif

extern MPI_Offset
ncmpii_hdr_len_NC(const NC *ncp, MPI_Offset sizeof_off_t);

extern int
ncmpii_hdr_get_NC(NC *ncp);

extern int 
ncmpii_hdr_put_NC(NC *ncp, void *buf);

extern int
ncmpii_NC_computeshapes(NC *ncp);

extern int
ncmpii_hdr_check_NC(bufferinfo *getbuf, NC *ncp);
/* end defined in header.c */

/* begin defined in mpincio.c */
extern int
ncmpiio_create(MPI_Comm comm, const char *path, int ioflags, MPI_Info info,
               ncio **nciopp);

extern int
ncmpiio_open(MPI_Comm comm, const char *path, int ioflags, MPI_Info info,
             ncio **nciopp);
extern int
ncmpiio_sync(ncio *nciop);

extern int
ncmpiio_move(ncio *const nciop, MPI_Offset to, MPI_Offset from,
             MPI_Offset nbytes);

extern int
ncmpiio_get_hint(NC *ncp, char *key, char *value, int *flag);

extern int
NC_computeshapes(NC *ncp);

/* end defined in mpincio.h */

/* begin defined in error.c */
const char * nc_strerror(int err);

void ncmpii_handle_error(int rank, int mpi_status, char *msg);
/* end defined in error.c */
/*
 * These functions are used to support
 * interface version 2 backward compatiblity.
 * N.B. these are tested in ../nc_test even though they are
 * not public. So, be careful to change the declarations in
 * ../nc_test/tests.h if you change these.
 */

extern int
ncmpii_put_att(int ncid, int varid, const char *name, nc_type datatype,
	       MPI_Offset len, const void *value);

extern int
ncmpii_get_att(int ncid, int varid, const char *name, void *value);

int ncmpii_x_putn_schar(void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_short(void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_int(void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_float(void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_putn_double(void *xbuf, const void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_schar(const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_short(const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_int(const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_float(const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);
int ncmpii_x_getn_double(const void *xbuf, void *buf, MPI_Offset nelems,
                        MPI_Datatype datatype);

int NCedgeck(const NC *ncp, const NC_var *varp, const MPI_Offset *start,
                const MPI_Offset *edges);

int NCstrideedgeck(const NC *ncp, const NC_var *varp, const MPI_Offset *start,
                const MPI_Offset *edges, const MPI_Offset *stride);

int NCcoordck(NC *ncp, const NC_var *varp, const MPI_Offset *coord);

int ncmpii_echar(nc_type nctype,MPI_Datatype mpitype);

int ncmpii_need_convert(nc_type nctype,MPI_Datatype mpitype);

int ncmpii_need_swap(nc_type nctype,MPI_Datatype mpitype);

void ncmpii_in_swapn(void *buf, MPI_Offset nelems, int esize);

int ncmpii_is_request_contiguous(NC_var *varp, const MPI_Offset starts[],
                const MPI_Offset  counts[]);

int ncmpii_get_offset(NC *ncp, NC_var *varp, const MPI_Offset starts[],
                const MPI_Offset counts[], const MPI_Offset strides[],
                MPI_Offset *offset_ptr);

int ncmpii_check_mpifh(NC* ncp, MPI_File *mpifh, MPI_Comm comm,
                int collective);

int ncmpii_update_numrecs(NC *ncp, MPI_Offset newnumrecs);

int ncmpii_vars_create_filetype(NC* ncp, NC_var* varp,
                const MPI_Offset start[], const MPI_Offset count[],
                const MPI_Offset stride[], int rw_flag,
                MPI_Offset *offset, MPI_Datatype *filetype);

extern int
ncmpii_getput_vars(NC *ncp, NC_var *varp, const MPI_Offset *start,
                const MPI_Offset *count, const MPI_Offset *stride,
                void *buf, MPI_Offset bufcount, MPI_Datatype datatype,
                int rw_flag, int io_method);

extern int
ncmpii_getput_varm(NC *ncp, NC_var *varp, const MPI_Offset start[],
		const MPI_Offset count[], const MPI_Offset stride[],
		const MPI_Offset imap[], void *buf, MPI_Offset bufcount,
		MPI_Datatype datatype, int rw_flag, int io_method);
extern int
ncmpii_igetput_varm(NC *ncp, NC_var *varp, const MPI_Offset *start,
                const MPI_Offset *stride, const MPI_Offset *imap,
                const MPI_Offset *count, void *buf, MPI_Offset bufcount,
                MPI_Datatype datatype, int *reqid, int rw_flag);

extern int
ncmpii_wait(NC *ncp, int io_method, int num_reqs, int *req_ids,
                int *statuses);

extern int
ncmpii_cancel(NC *ncp, int num_req, int *req_ids, int *statuses);

#endif /* _NC_H_ */
