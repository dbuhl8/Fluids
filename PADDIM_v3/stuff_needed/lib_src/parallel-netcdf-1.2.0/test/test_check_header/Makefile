srcdir		= .


INCDIR		= $(srcdir)/../../src/lib
LIBDIR		= ../../src/lib

INCLUDES	= -I$(INCDIR) -I../../src/lib -I$(srcdir)/../common

CC		= mpiicc
CFLAGS          = -O3 -xAVX $(INCLUDES) 

SRCS		= test_check_header.c test_check_header1.c
OBJS		= $(SRCS:.c=.o)
PROGS		= test_check_header test_check_header1
LDFLAGS		= -L$(LIBDIR) -L../common -xT
LIBS		= -lpnetcdf -ltestutils 

.c.o:
	$(CC) $(CFLAGS) -c $<

all:		$(PROGS)

test_check_header: test_check_header.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

test_check_header1: test_check_header1.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(PROGS) *.nc

distclean: clean
	rm -f Makefile
