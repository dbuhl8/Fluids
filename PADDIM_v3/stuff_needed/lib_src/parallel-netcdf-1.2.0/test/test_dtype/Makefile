srcdir		= .


INCDIR		= $(srcdir)/../../src/lib 
LIBDIR		= ../../src/lib

INCLUDES	= -I$(INCDIR) -I../../src/lib

CC		= mpiicc
CFLAGS          = -O3 -xAVX $(INCLUDES) 

SRCS		= test_subarray.c test_darray.c test_nonblocking.c
OBJS		= $(SRCS:.c=.o)
PROGS		= test_subarray test_darray test_nonblocking
LDFLAGS		= -L$(LIBDIR) 
LIBS		= -lpnetcdf 

.c.o:
	$(CC) $(CFLAGS) -c $<

all:		$(PROGS)

test_subarray: test_subarray.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

test_darray: test_darray.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

test_nonblocking: test_nonblocking.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(PROGS) *.nc

distclean: clean
	rm -f Makefile
