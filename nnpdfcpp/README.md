You need to ensure that:

 - libnnpdf
 - yaml-cpp
 - CERN-ROOT (for validphys)
 - APFEL

Are all working for you and can be found by your linker.

Then run:

cmake .

and then:

make

More options are available when running:

ccmake .

or

cmake-gui .