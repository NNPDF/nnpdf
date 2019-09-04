Introduction to `validphys 2`
=============================

`validphys 2` is a Python code that implements the data model of NNPDF
resources. It provides an executable, called `validphys` which is used to
analyze NNPDF specific data, which takes runcards written in
[YAML](https://en.wikipedia.org/wiki/YAML) as an input and can produce plots,
tables or entire reports as an output. The code also provides a Python library
(also called `validphys`) which is used to implement executables providing
interfaces to more specific analyses such as the `vp-comparefits`, and to
serve as basis to other NNPDF codes such as `n3fit`.

`validphys 2` is implemented on top of the `reportengine` framework.
`reportengine` provides the logic to process the runcards by building  task
execution graphs based on individual actions (which are Python functions). The
runcards can execute complex analysis and parameter scans with  the appropriate
use of namespaces.

Some parts of validphys use the `libnnpdf` library in C++, through SWIG
wrappers.

The ideas behind the design of the code are explained in the
[Design](./design.md) section.

