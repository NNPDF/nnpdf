default:
	mkdir -p pylibs/NNPDF
	cp src/*.py pylibs/NNPDF/
	find src -maxdepth 1 -name "*.i" -print0 | xargs -0 -I {} swig -c++ -outdir pylibs/NNPDF -python {}
	mkdir -p _sources
	mv src/*.cxx _sources

install:
	python setup.py install

clean:
	rm -rf _sources/
	rm -rf NNPDF/ build dist pylibs
