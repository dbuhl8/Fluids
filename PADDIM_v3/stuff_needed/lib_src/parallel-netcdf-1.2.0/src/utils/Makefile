srcdir = .

SUBDIRS = ncgen ncdump ncdiff

all: ncmpigen ncmpidump ncmpidiff

ncmpigen:
	cd ncgen && $(MAKE)

ncmpidump:
	cd ncdump && $(MAKE)

ncmpidiff:
	cd ncdiff && $(MAKE)

clean:
	set -e; for i in $(SUBDIRS); do ( cd $$i && $(MAKE) clean ; ) ; done

distclean:
	set -e; for i in $(SUBDIRS); do ( cd $$i && $(MAKE) distclean ; ) ; done
	rm -f Makefile
	if [ $(srcdir) != . ] ; then rmdir $(SUBDIRS) ; fi

install: ncgen_install ncdump_install ncdiff_install

ncgen_install:
	cd ncgen && $(MAKE) install

ncdump_install:
	cd ncdump && $(MAKE) install

ncdiff_install:
	cd ncdiff && $(MAKE) install
