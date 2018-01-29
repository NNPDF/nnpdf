################################################################################
#
# Makefile to compile and link C programs
#
# Version valid for Linux machines
#
# "make" compiles and links the specified main programs and modules
# using the specified libraries (if any), and produces the executables
# 
# "make clean" removes all files created by "make"
#
# Dependencies on included files are automatically taken care of
#
################################################################################

all: rmxeq mkdep mkxeq
	mv $(MAIN) ../bin
.PHONY: all


GCC = g++

# MACRO = -DHOPPET -DAPPLGRID

# main programs and required modules 

MAIN =  nnpdfrw

MODULES = pdfs utils

# search path for modules
#MDIR = .
#VPATH = .:$(MDIR)/PDF:$(MDIR)/appl

# root
ROOTINCS = $(shell root-config --cflags) 
ROOTLIBS = $(shell root-config --glibs) 

#LHAPDF
LHAPDFINCS = -I$(shell lhapdf-config --prefix)/include
LHAPDFDIR  = $(shell lhapdf-config --prefix)/lib
LHAPDFLIBS = -L$(LHAPDFDIR) -lLHAPDF

# additional include directories

#INCPATH = -I../../include $(ROOTINCS) $(LHAPDFINCS) 
INCPATH = $(ROOTINCS) $(LHAPDFINCS) 

# additional libraries to be included 

LIBS = $(ROOTLIBS) $(LHAPDFLIBS)

# scheduling and optimization options (such as -DSSE -DSSE2 -DP4)
 
CFLAGS = -Wall -ansi -O $(MACRO)

############################## do not change ###################################

SHELL=/bin/bash

CC=$(GCC)

PGMS= $(MAIN) $(MODULES)

INCDIRS = $(INCPATH)

OBJECTS = $(addsuffix .o,$(MODULES))

LDFLAGS = $(LIBS)

-include $(addsuffix .d,$(PGMS))


# rule to make dependencies

$(addsuffix .d,$(PGMS)): %.d: %.cc Makefile
	@ $(CC) -MM -ansi $(INCDIRS) $< -o $@


# rule to compile source programs

$(addsuffix .o,$(PGMS)): %.o: %.cc Makefile
	$(CC) $< -c $(CFLAGS) $(INCDIRS) -o $@


# rule to link object files

$(MAIN): %: %.o $(OBJECTS) Makefile
	$(CC) $< $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@


# produce executables

mkxeq: $(MAIN)


# remove old executables and old error log file

rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"		


# make dependencies

mkdep:  $(addsuffix .d,$(PGMS))
	@ echo "generate tables of dependencies"


# clean directory 

clean:
	@ -rm -rf *.d *.o *~ $(MAIN)

.PHONY: clean

################################################################################
