srcdir		= .


INCDIR		= $(srcdir)/../../src/lib -I../../src/lib
LIBDIR		= ../../src/lib

INCLUDES	= -I$(INCDIR) -I$(srcdir)/../common

CC		= mpiicc
CFLAGS          = -O3 -xAVX $(INCLUDES) 

SRCS		= test_inq_format.c 
OBJS		= $(SRCS:.c=.o) 
PROGS		= test_inq_format  
LDFLAGS		= -L$(LIBDIR) -L../common -xT
LIBS		= -lpnetcdf -ltestutils 

.c.o:
	$(CC) $(CFLAGS) -c $<

all:		$(PROGS)

test_inq_format: test_inq_format.o
	$(CC) -g -o $@ $@.o $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o $(PROGS) *.nc

distclean: clean
	rm -f Makefile
