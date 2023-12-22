srcdir		= .


INCDIR		= $(srcdir)/../../src/lib
LIBDIR		= ../../src/lib

INCLUDES	= -I$(INCDIR) -I$(srcdir)/../common -I../../src/lib

CC		= mpiicc
CFLAGS          = -O3 -xAVX $(INCLUDES) 

SRCS		= large_files.c
OBJS		= $(SRCS:.c=.o)
PROGS		= large_files
LDFLAGS		= -L$(LIBDIR) -L../common -xT
LIBS		= -lpnetcdf -ltestutils 

.c.o:
	$(CC) $(CFLAGS) -c $<

all:		$(PROGS)

large_files: large_files.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(PROGS) *.nc

distclean: clean
	rm -f Makefile
