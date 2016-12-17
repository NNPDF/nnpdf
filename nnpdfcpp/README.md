You need to ensure that:

 - libnnpdf
 - yaml-cpp
 - CERN-ROOT (for validphys)

Are all working for you and can be found by your linker.

Then run:

./autogen.sh

to generate the configure file, and then:

./configure
make
make install 
