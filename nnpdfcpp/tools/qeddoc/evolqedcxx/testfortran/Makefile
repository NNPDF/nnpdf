# Example of compilation with fortran

FF      = gfortran
FLAGS   = -O3
INCDIR  = -I../
LIBSDIR = -L../
LIBS    = -levolqed -lgsl -llapack -lstdc++

all: testevolqed.f
	make -C ../
	$(FF) $(FLAGS) -c testevolqed.f $(INCDIR)
	$(FF) $(FLAGS) -o testevolqed testevolqed.o $(LIBSDIR) $(LIBS)

clean:
	rm -rf *.o testevolqed