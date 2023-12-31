srcdir		= .


INCDIR		= $(srcdir)/../../src/lib 
LIBDIR		= ../../src/lib

INCLUDES	= -I$(INCDIR)  -I../../src/lib

CC		= mpiicc
CFLAGS          = -O3 -xAVX $(INCLUDES) 
FC		= mpiifort
FFLAGS          = -O3 -xAVX $(INCLUDES) -I../../src/libf
CXX		= 
CXXFLAGS        = -O3 -xAVX $(INCLUDES)


SRCS		= $(srcdir)/ncmpi_vars_null_stride.c $(srcdir)/vecors.c \
			$(srcdir)/bigrecords.F $(srcdir)/collective_error.c 
PROGS		= ncmpi_vars_null_stride vectors bigrecords redef1 collective_error
LDFLAGS		= -L$(LIBDIR) -xT
LIBS		= -lpnetcdf 

.c.o:
	$(CC) $(CFLAGS) -c $<
.F.o:
	$(FC) $(FFLAGS) -c $<

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

all:		$(PROGS)

ncmpi_vars_null_stride: ncmpi_vars_null_stride.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

vectors: vectors.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

bigrecords: bigrecords.o
	$(FC) -o $@ $@.o $(LDFLAGS) $(LIBS)

redef1: redef1.o
	$(CXX) -o $@ $@.o $(LDFLAGS) $(LIBS)

collective_error: collective_error.o
	$(CC) -o $@ $@.o $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o $(PROGS) *.nc

distclean: clean
	rm -f Makefile

test: ncmpi_vars_null_stride vectors
	./ncmpi_vars_null_stride test.nc && echo '*** ncmpi_vars_null_stride: Success ***'
	./vectors && echo '*** vectors: Success ***'

testing-cxx: redef1
	$(srcdir)/redef1.sh ./redef1  $(srcdir)/redef-good.ncdump && echo '*** redef1: Success ***'

testing-f: bigrecords
	./bigrecords && echo '*** bigrecords: Success ***'

testing-collective_error:
	echo "please run with 2 processes. e.g. mpiexec -n 2 collective_error"

testing: test
