EXTRA_PROGRAMS = gencovmat fkcheck mkthpredictions apfelevol tconvcheck

COMMON = ../common/src/nnpdfsettings.cc \
	 ../common/src/common.cc \
	 ../common/src/loadutils.cc \
	 ../common/src/datautils.cc \
         ../common/src/md5.cc \
	 ../common/src/theorydb.cc

gencovmat_SOURCES = src/gencovmat.cc $(COMMON)
gencovmat_CPPFLAGS = -Iinc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11
gencovmat_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(YAML_LDFLAGS)

fkcheck_SOURCES = src/fkcheck.cc $(COMMON)
fkcheck_CPPFLAGS = -Iinc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11
fkcheck_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(YAML_LDFLAGS)

mkthpredictions_SOURCES = src/mkthpredictions.cc $(COMMON)
mkthpredictions_CPPFLAGS = -Iinc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11
mkthpredictions_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(YAML_LDFLAGS)

apfelevol_SOURCES = src/apfelevol.cc $(COMMON)
apfelevol_CPPFLAGS = -Iinc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(APFEL_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11
apfelevol_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(APFEL_LDFLAGS) $(YAML_LDFLAGS)

tconvcheck_SOURCES = src/tconvcheck.cc $(COMMON)
tconvcheck_CPPFLAGS = -Iinc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11
tconvcheck_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(YAML_LDFLAGS)

test: gencovmat fkcheck mkthpredictions apfelevol tconvcheck
	cp gencovmat $(prefix)/bin
	cp fkcheck $(prefix)/bin
	cp mkthpredictions $(prefix)/bin
	cp apfelevol $(prefix)/bin
	cp tconvcheck $(prefix)/bin
testclean:
	rm -rf gencovmat fkcheck mkthpredictions apfelevol tconvcheck

.PHONY: test testclean
