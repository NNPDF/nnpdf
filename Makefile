#############################################
# Makefile to generate BUILDMASTER Program
# make => for main version
# make clean   => to clean project
#############################################

# Macros for results folders
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

CXXFLAGS=-Wall -g -I ./inc -I ./src -std=c++11 $(NNPDFCXX) $(GSLCXX) $(MACROS)
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

buildmaster_src = src/buildmaster.o \
                  src/buildmaster_utils.o \
	   	  		  src/common.o 

filter_src = $(wildcard filters/*.cc)
filter_obj = $(filter_src:.cc=.o)

buildmaster_src += $(filter_obj)

buildmaster: $(buildmaster_src)
	$(CXX) $(CFLAGS) $(buildmaster_src) -o buildmaster $(LDFLAGS)
