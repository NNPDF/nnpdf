## Dependencies

You need to ensure that:

 - libnnpdf
 - yaml-cpp
 - CERN-ROOT (for validphys)
 - APFEL

Are all working for you and can be found by your linker.
See version requirements in [conda-recipe/meta.yaml](https://github.com/NNPDF/nnpdfcpp/blob/master/conda-recipe/meta.yaml).

## Installation

Then run:
```bash
cmake .
make
```

More options are available when running:
```bash
ccmake .
```
or
```
cmake-gui .
```
