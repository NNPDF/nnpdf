SUBDIRS = projects
ACLOCAL_AMFLAGS = -I m4

test:
	cd testcode && $(MAKE) $(AM_MAKEFLAGS) test
testclean:
	cd testcode && $(MAKE) $(AM_MAKEFLAGS) testclean
.PHONY: test testclean

install-data-hook:
	rsync -r --exclude=.svn $(top_builddir)/config $(prefix)
	rsync -r --exclude=.svn $(top_builddir)/results $(prefix)
	rsync -r --exclude=.svn $(top_builddir)/cluster $(prefix)
	rsync -r --exclude=.svn $(top_builddir)/scripts/* $(prefix)/bin