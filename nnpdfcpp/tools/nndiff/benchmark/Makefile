CC=g++
CFLAGS=$(shell nnpdf-config --cppflags) $(shell lhapdf-config --cppflags)
LDFLAGS=$(shell nnpdf-config --ldflags) $(shell lhapdf-config --ldflags) -std=c++11

all: benchmark

benchmark: main.o nndiff.o
	g++ $(CFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cc
	g++ $(CFLAGS) -c $< -std=c++11

clean:
	rm -rf benchmark *.o *~
