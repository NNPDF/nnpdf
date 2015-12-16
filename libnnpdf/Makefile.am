lib_LTLIBRARIES = libnnpdf.la
libnnpdf_la_SOURCES = 	src/common.cc \
			src/commondata.cc \
			src/fkgenerator.cc \
			src/fastkernel.cc \
			src/fkset.cc \
			src/dataset.cc \
			src/experiments.cc \
			src/logger.cc \
			src/pdfset.cc \
			src/thpredictions.cc \
			src/parametrisation.cc \
			src/lhapdfset.cc \
			src/positivity.cc \
			src/randomgenerator.cc \
			src/utils.cc \
			src/nnmpi.cc

otherincludedir = $(includedir)/NNPDF
otherinclude_HEADERS =  src/NNPDF/common.h \
		        src/NNPDF/timer.h \
		        src/NNPDF/commondata.h \
		        src/NNPDF/experiments.h \
		        src/NNPDF/lhapdfset.h \
		        src/NNPDF/pdfset.h \
		        src/NNPDF/randomgenerator.h \
		        src/NNPDF/parametrisation.h \
		        src/NNPDF/utils.h \
		 	src/NNPDF/dataset.h \
			src/NNPDF/fkgenerator.h \
			src/NNPDF/fastkernel.h \
			src/NNPDF/fkset.h \
			src/NNPDF/logger.h \
			src/NNPDF/positivity.h \
			src/NNPDF/thpredictions.h \
                        src/NNPDF/nnmpi.h \
                        src/NNPDF/exceptions.h


libnnpdf_la_LDFLAGS = -version-info 0:0:0 $(GSL_LDFLAGS) $(SIMD_FLAGS) $(OPENMP_CFLAGS) $(LHAPDF_LDFLAGS)
libnnpdf_la_CPPFLAGS= $(GSL_CPPFLAGS) $(SIMD_FLAGS) $(OPENMP_CFLAGS) -fPIC $(LHAPDF_CPPFLAGS)


# Utility programs
bin_PROGRAMS = FKconvolute
FKconvolute_SOURCES = src/fkconvolute.cc
FKconvolute_CPPFLAGS = $(GSL_CPPFLAGS) $(SIMD_FLAGS) $(OPENMP_CFLAGS) -fPIC $(LHAPDF_CPPFLAGS)
FKconvolute_LDFLAGS = $(GSL_LDFLAGS) $(SIMD_FLAGS) $(OPENMP_CFLAGS) $(LHAPDF_LDFLAGS) -lnnpdf

# Scripts to install
dist_bin_SCRIPTS=scripts/FKmerge

ACLOCAL_AMFLAGS = -I m4
SUBDIRS = bin

wrapper:
	cd wrapper && $(MAKE) $(AM_MAKEFLAGS) install

wrapperclean:
	cd wrapper && $(MAKE) $(AM_MAKEFLAGS) clean

.PHONY: wrapper wrapperclean