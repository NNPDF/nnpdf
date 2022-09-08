.. _vp-index:

Code for data: validphys
========================

Introduction to ``validphys 2``
-------------------------------

* ``validphys 2`` is a Python code that implements the data model of NNPDF
  resources.

* It provides an executable, called ``validphys`` which is used to
  analyze NNPDF specific data, which takes runcards written in
  `YAML <https://en.wikipedia.org/wiki/YAML>`_ as an input and can produce plots,
  tables or entire reports as an output.

* The code also provides a Python library
  (also called `validphys`) which is used to implement executables providing
  interfaces to more specific analyses such as the ``vp-comparefits``, and to
  serve as basis to other NNPDF codes such as ``n3fit``.

* ``validphys 2`` is implemented on top of the
  `reportengine <https://github.com/NNPDF/reportengine>`_ framework.
  ``reportengine`` provides the logic to process the runcards by building task
  execution graphs based on individual actions (which are Python functions). The
  runcards can execute complex analysis and parameter scans with  the appropriate
  use of namespaces. More information on ``reportengine`` and its interface with
  ``validphys`` can be found in the
  :ref:`Design <design>` section.

* Some parts of ``validphys`` use the ``libnnpdf`` library in C++, through SWIG
  wrappers.

* The ideas behind the design of the code are explained in the
  :ref:`Design <design>` section.

Some things that `validphys` does
---------------------------------

* Download resources (``vp-get``) - see :ref:`download`
* Upload resources (``vp-upload``, ``wiki-upload`` and ``--upload`` flag) - see :ref:`upload`
* Prepare fits for running with ``n3fit`` (``vp_setupfit``) - see :ref:`scripts`
* Postprocess a fit (``postfit``) - see :ref:`scripts`
* Rename a fit or PDF (``vp-fitrename`` and ``vp-pdfrename``) - see :ref:`scripts`
* Sample a PDF (``vp-pdffromreplicas``) - see :ref:`scripts`
* Generate a report with information about possible inefficiencies in fitting methodology (``vp-deltachi2``) - see :ref:`scripts`
* Allow analysis via a high level interface - see :ref:`vpapi`
* Analyse results - see :ref:`tutorials`

Using validphys
---------------

.. toctree::
   :maxdepth: 1

   ./getting-started
   ./download
   ./upload
   ./nnprofile
   ./complex_runcards
   ./cuts
   ./datthcomp
   ./reports
   ./scripts

   ./developer
   ./tables_figs
   ./customplots
   ./examples

How validphys handles data
--------------------------
.. toctree::
   :maxdepth: 1

   ./pydataobjs
   ./filters
   ./theorycov/index
   ./dataspecification

More detailed functionality
---------------------------
.. toctree::
   :maxdepth: 1

   ./design
   ./namespaces
   ./resolving_dependencies
   ./checks
   ./collect
   ./api
   ./custom_pipelines
