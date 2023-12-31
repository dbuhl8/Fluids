# Makefile for netCDF (semi)exhaustive FORTRAN test.
#
# $Id: Makefile.in 816 2010-05-14 23:34:49Z wkliao $

srcdir = .


include ../../macros.make

# M4FLAGS	= -s -B7168
# CC		= cc -fullwarn -woff 1209,1506


INCLUDES	= -I$(srcdir) -I$(srcdir)/../../src/libf -I../../src/libf -I$(srcdir)/../../src/lib -I../../src/lib

FFLAGS	= $(INCLUDES) -O3 -xAVX 
lib_netcdf	= ../../src/lib/libpnetcdf.a
# Don't use "-L../libsrc -lnetcdf" in the following because that doesn't
# work on a CRAY T90 (sigh).
ld_netcdf	= $(lib_netcdf)

M4SRCS		= \
		  test_get.m4	\
		  test_put.m4   \
		  test_iget.m4  \
		  test_iput.m4
#

FSRCS		= \
		  $(M4SRCS:.m4=.F)	\
		  nf_error.F	\
		  nf_test.F	\
		  test_read.F	\
		  test_write.F	\
		  util.F
#

CSRCS		= fortlib.c

OBJS		= \
		  test_get.o	\
		  test_put.o	\
		  nf_error.o	\
		  nf_test.o	\
		  test_read.o	\
		  test_write.o	\
		  util.o	\
		  fortlib.o     \
		  test_iget.o   \
		  test_iput.o   

time_log	= times

GARBAGE		= \
		  $(time_log)	\
		  lint.out	\
		  nf_test	\
		  scratch.nc	\
		  *.f77		\
		  test.nc	\
		  $(M4SRCS:.m4=.F)
#

PACKING_LIST	= \
		  $(CSRCS)	\
		  $(FSRCS)	\
		  $(M4SRCS)	\
		  depend	\
		  Makefile	\
	          tests.inc
#

all:		nf_test

test:		nf_test test.nc
	./nf_test
	@echo '*** Success ***'

readonly:	nf_test test.nc
	./nf_test -r

test.nc:  nf_test
	./nf_test -c

install:

uninstall:


nf_test:	$(OBJS) $(lib_netcdf)
	$(LINK.F) $(OBJS) $(ld_netcdf) $(FLIBS) $(MATHLIB) $(LIBS)


tags:		$(FSRCS) $(CSRCS) FORCE
	ctags -t $(FSRCS) $(CSRCS) ../fortran/*.c ../libsrc/*.c


test_get.F:	test_get.m4
test_put.F:	test_put.m4
test_iget.F:	test_iget.m4
test_iput.F:	test_iput.m4

# This simple testing target ensures that the test files are present
testing: test.nc nf_test
	rm -f scratch.nc
	./nf_test

#test.nc: 
#	(cd ../nc_test && make nc_test && ./nc_test -c)
#	cp ../nc_test/test.nc .

include $(srcdir)/../../rules.make
include $(srcdir)/depend
