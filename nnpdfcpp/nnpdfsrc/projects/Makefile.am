bin_PROGRAMS = filter chi2check validphys nnfit

COMMON = ../common/src/common.cc \
	 ../common/src/nnpdfsettings.cc \
	 ../common/src/loadutils.cc \
         ../common/src/datautils.cc \
      ../common/src/md5.cc \
         ../common/src/theorydb.cc

filter_SOURCES = filtersrc/src/filter.cc \
	         filtersrc/src/kincuts.cc $(COMMON)
filter_CPPFLAGS = -Ifiltersrc/inc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11
filter_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(YAML_LDFLAGS)

chi2check_SOURCES = chi2checksrc/src/chi2check.cc $(COMMON)
chi2check_CPPFLAGS = -Ichi2checksrc/inc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11
chi2check_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(YAML_LDFLAGS)

nnfit_SOURCES = nnfitsrc/src/nnfit.cc \
		nnfitsrc/src/fitbases.cc \
		nnfitsrc/src/pdfbasis.cc \
		nnfitsrc/src/minimizer.cc \
		nnfitsrc/src/fastaddchi2.cc \
		nnfitsrc/src/fitpdfset.cc \
		nnfitsrc/src/stopping.cc \
		nnfitsrc/src/apfelevol.cc \
		$(COMMON)
nnfit_CPPFLAGS = -Innfitsrc/inc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11 $(APFEL_CPPFLAGS)
nnfit_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(YAML_LDFLAGS) $(APFEL_LDFLAGS)

validphys_SOURCES = validphyssrc/src/validphys.cc \
		    validphyssrc/src/plotdata.cc \
		    validphyssrc/src/plotutils.cc $(COMMON)
validphys_CPPFLAGS = -Ivalidphyssrc/inc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11 $(ROOT_CPPFLAGS)
validphys_LDFLAGS = $(NNPDF_LDFLAGS) $(ROOT_LDFLAGS) $(SQLITE3_LDFLAGS) $(YAML_LDFLAGS)
