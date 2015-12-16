
if WITH_PYWRAP

PYEXT_ENV = CC=$(CC) CXX=$(CXX)

all-local:
	mkdir -p pylibs/NNPDF
	cp src/*.py pylibs/NNPDF/
	find src -name "*.i" -print0 | xargs -0 -I {} $(SWIG) -c++ -outdir pylibs/NNPDF -python {}
	mkdir -p _sources
	mv src/*.cxx _sources
	$(PYEXT_ENV) $(PYTHON) setup.py install

install-exec-local:
	$(PYEXT_ENV) $(PYTHON) setup.py install --force

clean-local:
	-rm -rf $(builddir)/_sources/
	-rm -rf $(builddir)/NNPDF/ 
	-rm -rf $(builddir)/build 
	-rm -rf $(builddir)/dist 
	-rm -rf $(builddir)/pylibs

distclean-local: clean-local

endif
