# NNPDF flags
NFLAGS=$(shell nnpdf-config --cppflags)
NLDFLAGS=$(shell nnpdf-config --ldflags)

# LHAPDF flags
LHAFLAGS= $(shell lhapdf-config --cppflags)
LHALDFLAGS=  $(shell lhapdf-config --ldflags)


#########################################

# Total flags for test code
TSTFLAGS=$(NFLAGS) $(LHAFLAGS)
TSTLDFLAGS=$(NLDFLAGS) $(LHALDFLAGS)

#########################################

.PHONY: clean install

all:	test

test:
	g++ -Wall $(TSTFLAGS) $(TSTLDFLAGS) fkconvolute.cc -o fkconvolute
	
clean:
	rm -f fkconvolute
