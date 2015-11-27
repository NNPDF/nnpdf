#############################################
# Makefile to generate BUILDMASTER Program
# make release => for release version
# make debug   => for debug version
# make doc     => for doxygen documentation
# make clean   => for clean project
#############################################

SOURCEDIR=./src
INCLUDEDIR=inc
OBJECTDIR=obj
DESTDIR=.

VPATH=$(SOURCEDIR)
vpath %.h $(COREINC)
vpath %.cc $(CORESRC)

# refers to bin folder
RESULTDIR = results
RESULTSDIR= -D  RESULTS_PATH="results"
DATADIR=    -D  DATA_PATH="./"

# Compiler options
##################

CC=g++
CFLAGS=-Wall -O2 -Wno-unknown-pragmas
LDFLAGS=-c -Wall -O2 -Wno-unknown-pragmas
debug: CFLAGS=-Wall -O0 -g -Wno-unknown-pragmas
debug: LDFLAGS=-c -Wall -O0 -g -Wno-unknown-pragmas

GSLINCLUDE= $(shell gsl-config --cflags)
GSLLIBS= $(shell gsl-config --libs)

NNPDFINCLUDE = $(shell nnpdf-config --cppflags) -I$(INCLUDEDIR) -std=c++11
NNPDFLIBS = $(shell nnpdf-config --ldflags)
NNPDFMACROS = $(DATADIR)

#####################################

all: release

release: svn $(DESTDIR)/buildmaster
debug: svn $(DESTDIR)/buildmaster


svn:
	rm -rf $(INCLUDEDIR)/svn.h; \
	echo "#define SVN_REV \"$(shell svnversion -n)\"" > $(INCLUDEDIR)/svn.h

uninstall: clean

clean:
	rm -rf $(DESTDIR)/buildmaster
	rm -rf $(OBJECTDIR)/*.*
	rm -rf $(RESULTDIR)/*

######### Programs BUILDMASTER #########
########################################

buildmaster_obj = buildmaster.o \
                  buildmaster_utils.o \
	   	  common.o \
		  NMC.o \
		  SLAC.o  \
	    	  BCDMS.o \
                  ATLAS.o \
                  ATLAS2011JETS.o \
                  ATLASLOMASSDY11.o \
                  CMS.o   \
                  CDF.o   \
                  LHCb.o  \
                  D0.o    \
                  FTDY.o  \
		  		  H1HERA2.o \
		  		  HERACOMB.o \
                  CHORUS.o \
                  NUTEV.o \
                  HERA1-C.o \
                  HERA2-C.o \
                  ZEUS2.o \
                  TOP.o \
                  CMSwc.o \
                  CMSDY2D12.o \
                  ATLASTOPDIFF.o \
				  CMSTOPDIFF.o \
                    POS.o \

buildmaster_src = $(addprefix $(OBJECTDIR)/, $(buildmaster_obj))

$(DESTDIR)/buildmaster: $(buildmaster_src)
	$(CC) $(CFLAGS) $(NNPDFMACROS) $(NNPDFINCLUDE) $^ $(NNPDFLIBS) \
	-o $@
	@echo "==> $@ done!"

## o libraries

$(OBJECTDIR)/%.o: %.cc
	$(CC) $(NNPDFMACROS) $(LDFLAGS) $(NNPDFINCLUDE) $< -o $@
