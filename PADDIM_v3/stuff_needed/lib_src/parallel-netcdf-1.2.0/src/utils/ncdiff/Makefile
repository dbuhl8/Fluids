# Makefile for ncmpidump(1).
#
# $Id: Makefile.in 762 2010-01-20 22:55:07Z robl $

srcdir = .

INSTALL = /usr/bin/install -c

include ../../../macros.make

NCGEN		= ../ncgen/ncmpigen
INCLUDES	= -I$(srcdir)/../../../src/lib -I../../lib -I$(srcdir) 

c_sources	= ncmpidiff.c
headers		=

PROGRAM		= ncmpidiff

PACKING_LIST	= $(c_sources) $(headers) depend \
		  ncmpidiff.1 Makefile
MANUAL		= ncmpidiff.1

lib_netcdf	= ../../../src/lib/libpnetcdf.a
ld_netcdf	= -L../../../src/lib -lpnetcdf

OBJS		=  ncmpidiff.o
GARBAGE		= $(PROGRAM)

all:		$(PROGRAM)

$(PROGRAM):	$(lib_netcdf) $(OBJS)
	$(LINK.c) $(OBJS) $(ld_netcdf) $(LIBS) 

install:	
	$(INSTALL) -d -m 755 $(MANDIR)/manl
	$(INSTALL) -m 644 $(srcdir)/$(MANUAL) $(MANDIR)/man1/$(MANUAL)

	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 755 $(PROGRAM) $(BINDIR)/$(PROGRAM)

uninstall:
	-rm -f $(BINDIR)/$(PROGRAM)
	-rm -f $(MANDIR)/man1/$(MANUAL)

$(PROGRAM)_oc : $(c_sources)
	#setopt primary_language C
	#load -C $(CPPFLAGS) $(c_sources)
	#load -C $(LIBS)
	#setopt program_name $(PROGRAM)

TAGS:		FORCE
	etags `echo $(PACKING_LIST) | fmt -1 | egrep '\.c|\.h'

include $(srcdir)/../../../rules.make
include $(srcdir)/depend
