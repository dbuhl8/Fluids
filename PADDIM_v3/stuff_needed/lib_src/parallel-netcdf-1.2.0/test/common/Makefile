srcdir = .


CC = mpiicc
CFLAGS = -O3 -xAVX

OBJS = testutils.o

all: libtestutils.a

libtestutils.a: $(OBJS)
	ar cru libtestutils.a $(OBJS)
testutils.o: testutils.c testutils.h

clean:
	rm -f *.o *.a

distclean: clean
	rm -f Makefile
