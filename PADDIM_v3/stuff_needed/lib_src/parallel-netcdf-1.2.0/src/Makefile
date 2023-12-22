srcdir = .


all:
	cd lib && $(MAKE)
	if [ yes = yes ] ; then cd libf && $(MAKE) ; fi
	cd utils && $(MAKE)

install:
	cd lib && $(MAKE) install
	if [ yes = yes ] ; then cd libf && $(MAKE) install ; fi
	cd utils && $(MAKE) install

clean:
	cd lib && $(MAKE) clean
	if [ yes = yes ] ; then cd libf && $(MAKE) clean ; fi
	cd utils && $(MAKE) clean

distclean:
	cd lib && $(MAKE) distclean
	if [ yes = yes ] ; then cd libf && $(MAKE) distclean ; fi
	cd utils && $(MAKE) distclean
	rm -f Makefile
	if [ $(srcdir) != . ] ; then if [ yes = yes ] ; then rmdir libf ; fi ; fi
	if [ $(srcdir) != . ] ; then rmdir lib utils ; fi
