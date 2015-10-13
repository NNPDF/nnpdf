#############################################
# Makefile to generate c++ programs
# make release => for release version
# make debug   => for debug version
# make doc     => for doxygen documentation
# make clean   => for clean project
#############################################

# Files
##################
EXECUTABLE= comparison

EXECUTABLE_DEBUG= 

CCFILES = 

# Compiler options
##################
CC=g++
CFLAGS=-Wall -O3 $(shell root-config --cflags)
CFLAGSD=-Wall -O0 -g $(shell root-config --cflags)
LDFLAGS=-c -Wall -O3 $(shell root-config --cflags)
LDFLAGSD=-c -Wall -O0 -g $(shell root-config --cflags)
LIBS= -L$(shell lhapdf-config --libdir) -lLHAPDF $(shell root-config --glibs)
INCL= -I$(shell lhapdf-config --incdir)

# Data
##################
SRCDIR=
INCLUDEDIR=
FILENAME=$(addsuffix .cc, $(CCFILES))
FILENAMEH=$(addsuffix .h, $(HFILES))
SOURCES=$(addprefix $(SRCDIR), $(FILENAME))
INCLUDES=$(addprefix $(INCLUDEDIR), $(FILENAMEH))
OBJECTS=$(SOURCES:.cc=.o)
OBJECTSD=$(SOURCES:.cc=_d.o)

FILEMAIN=$(addsuffix .cc, $(EXECUTABLE))
SOURCEMAIN=$(addprefix $(SRCDIR), $(FILEMAIN))
MAINS=$(SOURCEMAIN:.cc=.o)
MAINSD=$(SOURCEMAIN:.cc=_d.o)

all: release

install: all

uninstall: clean

release: $(SOURCES) $(MAINS) $(EXECUTABLE)

debug: $(SOURCES) $(MAINSD) $(EXECUTABLE_DEBUG)

$(EXECUTABLE_DEBUG): $(OBJECTSD) $(INCLUDES) $(SRCDIR)/$<
	$(CC) $(CFLAGSD) $(OBJECTSD) $(SRCDIR)$@.o -I $(INCLUDEDIR) $(LIBS) $(INCL) \
	-o $@
	@echo "==> Build $@"

%_d.o: %.cc 
	$(CC) $(LDFLAGSD) -I $(INCLUDEDIR) $(INCL) $(LIBS) $< -o $@

$(EXECUTABLE): $(OBJECTS) $(INCLUDES) $(SRCDIR)/$<
	$(CC) $(CFLAGS) $(OBJECTS) $(SRCDIR)$@.o -I $(INCLUDEDIR) $(INCL) $(LIBS) \
	-o $@
	@echo "==> Build $@"

%.o: %.cc
	$(CC) $(LDFLAGS) -I $(INCLUDEDIR) $(LIBS) $(INCL) $< -o $@

clean:
	rm -rf $(EXECUTABLE) $(EXECUTABLE_DEBUG)
	rm -rf $(SRCDIR)*.o $(SRCDIR)*.d $(SRCDIR)*~ *.eps *.root  *.C *.rep* comparison 
