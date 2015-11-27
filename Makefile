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

# refers to bin folder
RESULTDIR = results
RESULTSDIR= -D  RESULTS_PATH="results"
DATADIR=    -D  DATA_PATH="./"
MACROS = $(RESULTSDIR) $(DATADIR)

# Compiler options
##################

CXX=c++

# NNPDF flags
NNPDFCXX = $(shell nnpdf-config --cppflags)
NNPDFLD = $(shell nnpdf-config --ldflags)

# GSL flags
GSLCXX = $(shell gsl-config --cflags)
GSLLD = $(shell gsl-config --libs)

CXXFLAGS=-Wall -I ./inc -std=c++11 $(NNPDFCXX) $(GSLCXX) $(MACROS)
LDFLAGS= $(NNPDFLD) $(GSLLD)

#####################################

all: buildmaster

.PHONY: clean
clean:
	rm -rf buildmaster
	rm -rf ./src/*.o
	rm -rf ./filters/*.o
	rm -rf $(RESULTDIR)/*

######### Programs BUILDMASTER #########
########################################

buildmaster_src = $(SOURCEDIR)/buildmaster.o \
                  $(SOURCEDIR)/buildmaster_utils.o \
	   	  		  $(SOURCEDIR)/common.o 

filter_src = $(wildcard filters/*.cc)
filter_obj = $(filter_src:.cc=.o)

buildmaster_src += $(filter_obj)

buildmaster: $(buildmaster_src)
	$(CXX) $(CFLAGS) $(buildmaster_src) -o buildmaster $(LDFLAGS)
