# To run v6
# export DYLD_LIBRARY_PATH=/Users/juanrojo/software/lhapdf/lha-6.1.3/lib/
# export LD_LIBRARY_PATH=/Users/juanrojo/software/lhapdf/lha-6.1.3/lib/


lhapdftest_v6:
	g++  -I/Users/juanrojo/software/lhapdf/lha-6.1.3/include  -L/Users/juanrojo/software/lhapdf/lha-6.1.3/lib -lLHAPDF   lhapdftest_v6.cc -o lhapdftest_v6.o
	cp lhapdftest_v6.o lhapdftest_v6

# to run v6
# export DYLD_LIBRARY_PATH=/Users/juanrojo/software/lhapdf/lha-5.9.0/lib/
# export LD_LIBRARY_PATH=/Users/juanrojo/software/lhapdf/lha-5.9.0/lib/

lhapdftest_v5:
	g++  -I/Users/juanrojo/software/lhapdf/lha-5.9.0/include  -L/Users/juanrojo/software/lhapdf/lha-5.9.0/lib -lLHAPDF   lhapdftest_v5.cc -o lhapdftest_v5.o
	cp lhapdftest_v5.o lhapdftest_v5



clean:
	rm -rf  lhapdftest_v6.o lhapdftest_v6 *.dat *~ lhapdftest_v5.o lhapdftest_v5 *.res

