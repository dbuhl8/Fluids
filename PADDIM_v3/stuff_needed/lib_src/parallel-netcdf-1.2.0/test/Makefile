srcdir = .

CXX = @CXX@
SUBDIRS = common test_double test_double_int test_float test_int test_dtype nc_test fandc nf_test largefile C mcoll_perf test_inq_format test_check_header testcases

all:
	set -e; for i in $(SUBDIRS); do ( cd $$i && $(MAKE) ; ) ; done

clean:
	set -e; for i in $(SUBDIRS); do ( cd $$i && $(MAKE) clean ; ) ; done

distclean:
	set -e; for i in $(SUBDIRS); do ( cd $$i && $(MAKE) distclean ; ) ; done
	rm -f Makefile
	if [ . != . ] ; then rmdir $(SUBDIRS) ; fi

# Eventually, this target will run all of the relevant tests
testing:
	cd common && make
	cd C && make testing
	cd mcoll_perf && make testing
	cd nc_test && make testing
	if [ yes = yes ] ; then cd nf_test && make testing ; fi
	if [ yes = yes ] ; then cd mcoll_perf && make testing-f ; fi
	if [ yes = yes ] ; then cd testcases && make testing-f ; fi
	if [ $(CXX) != "@CXX@" ] ; then cd testcases && make testing-cxx ; fi

