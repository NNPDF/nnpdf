EXTRA_PROGRAMS = gencovmat fkcheck mkthpredictions tconvcheck revolve

COMMON = ../common/src/nnpdfsettings.cc \
	 ../common/src/common.cc \
	 ../common/src/loadutils.cc \
	 ../common/src/datautils.cc \
         ../common/src/md5.cc

gencovmat_SOURCES = src/gencovmat.cc $(COMMON)
gencovmat_CPPFLAGS = -Iinc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11
gencovmat_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(YAML_LDFLAGS)

fkcheck_SOURCES = src/fkcheck.cc $(COMMON)
fkcheck_CPPFLAGS = -Iinc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11
fkcheck_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(YAML_LDFLAGS)

mkthpredictions_SOURCES = src/mkthpredictions.cc $(COMMON)
mkthpredictions_CPPFLAGS = -Iinc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11
mkthpredictions_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(YAML_LDFLAGS)

tconvcheck_SOURCES = src/tconvcheck.cc $(COMMON)
tconvcheck_CPPFLAGS = -Iinc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11
tconvcheck_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(YAML_LDFLAGS)

revolve_SOURCES = src/revolve.cc $(COMMON)
revolve_CPPFLAGS = -Iinc $(NNPDF_CPPFLAGS) -I$(top_builddir)/common/inc $(LHAPDF_CPPFLAGS) $(SQLITE3_CPPFLAGS) $(APFEL_CPPFLAGS) $(YAML_CPPFLAGS) -std=c++11
revolve_LDFLAGS = $(NNPDF_LDFLAGS) $(SQLITE3_LDFLAGS) $(APFEL_LDFLAGS) $(YAML_LDFLAGS)

test: gencovmat fkcheck mkthpredictions tconvcheck revolve
	cp gencovmat $(prefix)/bin
	cp fkcheck $(prefix)/bin
	cp mkthpredictions $(prefix)/bin
	cp tconvcheck $(prefix)/bin
	cp revolve $(prefix)/bin

testclean:
	rm -rf gencovmat fkcheck mkthpredictions tconvcheck revolve

.PHONY: test testclean
