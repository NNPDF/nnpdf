
if WITH_PYWRAP

PYEXT_ENV = CC=$(CC) CXX=$(CXX)

all:
	mkdir -p pylibs/NNPDF
	cp src/*.py pylibs/NNPDF/
	find src -maxdepth 1 -name "*.i" -print0 | xargs -0 -I {} $(SWIG) -c++ -outdir pylibs/NNPDF -python -py3 {}
	mkdir -p _sources
	mv src/*.cxx _sources

install: all
	$(PYEXT_ENV) $(PYTHON) setup.py install

clean-local:
	-rm -rf $(builddir)/_sources/
	-rm -rf $(builddir)/NNPDF/ 
	-rm -rf $(builddir)/build 
	-rm -rf $(builddir)/dist 
	-rm -rf $(builddir)/pylibs

distclean-local: clean-local

endif
