srcdir = .


include ../../macros.make

ALL: pnctestf pnctest csnap pnf_test

INCDIR = ../../src/lib
INCDIRF = ../../src/libf/
LNKDIR = ../../src/lib

MPICC = mpiicc
MPIF77 = mpiifort

EXECS = pnctestf pnctest csnap pnf_test

pnctestf: pnctestf.F
	$(LINK.F) -o pnctestf $(srcdir)/pnctestf.F -I$(INCDIRF) -L$(LNKDIR) -lpnetcdf -lm $(LIBS)

pnctest: pnctest.c
	$(LINK.c) -o pnctest $(srcdir)/pnctest.c -I$(INCDIR) -L$(LNKDIR) -lpnetcdf -lm $(LIBS)

csnap:  csnap.o
	$(LINK.c) -L$(LNKDIR) -o csnap csnap.o -lpnetcdf -lm $(LIBS)

csnap.o:  csnap.c
	$(MPICC) -I$(INCDIR) -c $(srcdir)/csnap.c

pnf_test:  pnf_test.o
	$(LINK.F) -L$(LNKDIR) -o pnf_test pnf_test.o -lpnetcdf -lm $(LIBS)

pnf_test.o:  pnf_test.F
	$(MPIF77) $(FFLAGS) -I$(INCDIRF) -c $(srcdir)/pnf_test.F 

clean:
	rm -f $(EXECS) *.o *.nc

distclean: clean
	rm -f Makefile

