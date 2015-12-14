default:
	mkdir -p pylibs/NNPDF
	touch pylibs/NNPDF/__init__.py
	find src -name "*.i" -print0 | xargs -0 -I {} swig -c++ -outdir pylibs/NNPDF -python {}
	mkdir -p _sources
	mv src/*.cxx _sources

clean:
	rm _sources/*.cxx
	rm -r NNPDF/*
