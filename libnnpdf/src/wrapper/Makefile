default:
	mkdir -p NNPDF
	touch NNPDF/__init__.py
	find . -name "*.i" -print0 | xargs -0 -I {} swig -c++ -outdir NNPDF -python {}
	mkdir -p _sources
	mv *.cxx _sources

clean:
	rm _sources/*.cxx
	rm -r NNPDF/*
