# APPLgrid flags
AFLAGS=$(shell applgrid-config --cxxflags)
ALDFLAGS=$(shell applgrid-config --ldflags)

# LHAPDF flags
LHAFLAGS= $(shell lhapdf-config --cppflags)
LHALDFLAGS=  $(shell lhapdf-config --ldflags)


#########################################

# Total flags for test code
TSTFLAGS=$(AFLAGS) $(LHAFLAGS)
TSTLDFLAGS=$(ALDFLAGS) $(LHALDFLAGS)

#########################################

.PHONY: clean install

all:	test

test:
	g++ -Wall $(TSTFLAGS) $(TSTLDFLAGS) applgrid-test-pp.cpp -o applgrid-test-pp
	g++ -Wall $(TSTFLAGS) $(TSTLDFLAGS) applgrid-test-ppbar.cpp -o applgrid-test-ppbar

clean:
	rm -f applgrid-test-pp
	rm -f applgrid-test-ppbar
