srcdir		= .


INCDIR		= $(srcdir)/../../src/lib
LIBDIR		= ../../src/lib

INCLUDES	= -I$(INCDIR) -I$(srcdir)/../common -I../../src/lib

CC		= mpiicc
CFLAGS          = -O3 -xAVX $(INCLUDES) 

SRCS		= test_read.c test_read_indep.c test_write.c test_write_indep.c test_write64.c test_read64.c
OBJS		= $(SRCS:.c=.o)
PROGS		= test_read test_read_indep test_write test_write_indep test_write64 test_read64 
LDFLAGS		= -L$(LIBDIR) -L..//common -xT
LIBS		= -lpnetcdf -ltestutils 

.c.o:
	$(CC) $(CFLAGS) -c $<

all:		$(PROGS)

test_read: test_read.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)
test_read_indep: test_read_indep.o 
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)
test_write: test_write.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)
test_write_indep: test_write_indep.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)
test_write64: test_write64.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)
test_read64: test_read64.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(PROGS) *.nc

distclean: clean
	rm -f Makefile
