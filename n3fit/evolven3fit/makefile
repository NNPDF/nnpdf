## Simple and stupid makefile

NNPDFSRC=../../
INCLUDES=`pkg-config nnpdf --cflags` `pkg-config gsl --cflags` `apfel-config --cxxflags` `pkg-config yaml-cpp --cflags`
LIBS=`pkg-config nnpdf --libs` `pkg-config gsl --libs` `apfel-config --ldflags` `pkg-config yaml-cpp --libs`

all: evolven3fit

evolven3fit: evolven3fit.cc
	$(CXX) -I./ -I$(NNPDFSRC)/nnpdfcpp/src/nnfit/inc \
	$(INCLUDES) -I$(NNPDFSRC)/nnpdfcpp/src/common/inc \
	$(NNPDFSRC)/nnpdfcpp/src/common/src/nnpdfsettings.cc \
	$(NNPDFSRC)/nnpdfcpp/src/common/src/md5.cc \
	$(NNPDFSRC)/nnpdfcpp/src/common/src/exportgrid.cc \
	$(NNPDFSRC)/nnpdfcpp/src/nnfit/src/evolgrid.cc \
	$< -o $@ $(LIBS)

clean:
	rm -rf evolven3fit
