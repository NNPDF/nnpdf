
CXX = g++
F77 = gfortran

FFLAGS   += -O3 -fPIC 
CXXFLAGS += -O3 -fPIC 

# root
ROOTINCS = -I/$(shell root-config --cflags) 
ROOTLIBS = $(shell root-config --glibs) 
ROOTARCH = $(findstring -m64, $(ROOTINCS) )

# GSL
GSLINCS = $(shell gsl-config --prefix)/include 
GSLLIBS = $(shell gsl-config --libs) 

#LHAPDF
LHAPDFINCS = -I$(shell lhapdf-config --prefix)/include
LHAPDFDIR  = $(shell lhapdf-config --prefix)/lib
LHAPDFLIBS = -L$(LHAPDFDIR) -lLHAPDF

# applgrid
APPLCXXFLAGS = $(shell applgrid-config --cxxflags)
APPLCLIBS    = $(shell applgrid-config --ldcflags)
APPLFLIBS    = $(shell applgrid-config --ldflags)

# hoppet
HOPPETLIBS =  $(shell hoppet-config --libs)

# get the fotran runtime library for linking fortran 

FRTLLIB = $(shell gfortran $(CXXFLAGS) -print-file-name=libgfortran.a)
FRTLIB  = -L$(subst /libgfortran.a, ,$(FRTLLIB) ) -lgfortran

# now set up the compile and link flags and libs
CXXFLAGS += $(ROOTARCH) $(ROOTINCS) $(APPLCXXFLAGS) $(LHAPDFINCS) $(GSLINCS) 
LDFLAGS  += $(ROOTARCH)
FFLAGS   += $(ROOTARCH)

CLIBS += $(ROOTLIBS) $(LHAPDFLIBS) $(HOPPETLIBS) $(APPLCLIBS) $(GSLLIBS)
FLIBS += $(ROOTLIBS) $(LHAPDFLIBS) $(HOPPETLIBS) $(APPLFLIBS) $(FRTLIB)

install : all
all : distances

distances: distances.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS) 


.SUFFIXES : .cxx .o .f .c

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 

clean:
	rm -rf ./.libs ./.obj *.lo *.o *.la distances *.eps *~ *.tgz *.res

