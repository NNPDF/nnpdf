.. _reproduce40:

How to reproduce an NNPDF4.0 fit
================================================================================

Here we describe how to reproduce the NNPDF4.0 PDF sets that are publicly
available through LHAPDF as discussed in section 10 of the NNPDF4.0 paper
:cite:p:`nnpdf40`.


The files needed to reproduce the results are available in the folder
|n3fit_nnpdf40_folder|_
of the `project repository <https://github.com/NNPDF/nnpdf>`_ on Github. The
`.yml` files that folder are the following:

.. _nnpdf40 runcard textblock:

.. code-block::

    nnpdf40_env.yml # conda environment

    NNPDF40_hyperopt.yml # hyperoptimization runcard

    # Fit runcards:
    NNPDF40_lo_as_0118_pch.yml
    NNPDF40_lo_as_0118.yml
    NNPDF40_nlo_as_0117.yml
    NNPDF40_nlo_as_0118_pch.yml
    NNPDF40_nlo_as_0118.yml
    NNPDF40_nlo_as_0119.yml
    NNPDF40_nnlo_as_0116.yml
    NNPDF40_nnlo_as_01175.yml
    NNPDF40_nnlo_as_0117.yml
    NNPDF40_nnlo_as_0118_1000.yml
    NNPDF40_nnlo_as_01185.yml
    NNPDF40_nnlo_as_0118_pch.yml
    NNPDF40_nnlo_as_0119.yml
    NNPDF40_nnlo_as_0120.yml

    NNPDF40_nnlo_pdfas.yml
    NNPDF40_nnlo_as_0118_hessian.yml # Hessian conversion runcard


.. |n3fit_nnpdf40_folder| replace:: ``n3fit/runcards/reproduce_nnpdf40``
.. _n3fit_nnpdf40_folder: https://github.com/NNPDF/nnpdf/tree/master/n3fit/runcards/reproduce_nnpdf40

Below we will describe how these runcards can be used to generate the NNPDF4.0
PDF grids.

Setting up the NNPDF4.0 conda envirnoment
--------------------------------------------------------------------------------

The exact :ref:`conda environment <conda>`, including all transitive
dependencies, used to produce all the publicly released PDF fits has been
preserved to ensure precise reproducibility.  The environment can be found in
the `project repository
<https://github.com/NNPDF/nnpdf/tree/master/n3fit/runcards/reproduce_nnpdf40/nnpdf40_env.yml>`_
It can be generated and activated by running

.. code:: bash

    conda env create --file nnpdf40_env.yaml
    conda activate nnpdf40

Alternatively, we provide a ready to run docker image with a pre-installed
environment. After :ref:`setting it up <docker>`, you can download and open
the environment by running

.. code:: bash

    docker run -it ghcr.io/nnpdf/nnpdf:4.0.3 bash



Hyperoptimization for NNPDF4.0
--------------------------------------------------------------------------------
The |NNPDF40_hyperopt.yml|_ runcard can be used to run the
:ref:`hyperoptimization <hyperoptimization>` with the exact same settings used
to identify the NNPDF4.0 architecture. This step can be skipped if the aim is
simply to reproduce the NNPDF4.0 grids, which requires using the exact same
architecture as was used to create the NNPDF4.0 grids. Nevertheless, this
runcard is made available to provide a complete overview of each step in the
production of the NNPDF4.0 sets.

.. |NNPDF40_hyperopt.yml| replace:: ``NNPDF40_hyperopt.yml``
.. _NNPDF40_hyperopt.yml: https://github.com/NNPDF/nnpdf/tree/master/n3fit/runcards/reproduce_nnpdf40/NNPDF40_hyperopt.yml


Producing the NNPDF4.0 fits
--------------------------------------------------------------------------------
After setting up the conda environment, we can use it to produce the fits in the
:ref:`usual way <run-n3fit-fit>`. The names of the "fit runcards" listed in the
:ref:`textblock <nnpdf40 runcard textblock>` above are the same as those of
the corresponding (public) LHAPDF grids of the same name. For example,
baseline fit with :math:`\alpha_s(m_Z)=0.118` and a  variable-flavor-number
scheme with up to five active flavors can be generated with
|NNPDF40_nnlo_as_0118_1000.yml|_.

.. |NNPDF40_nnlo_as_0118_1000.yml| replace:: ``NNPDF40_nnlo_as_0118_1000.yml``
.. _NNPDF40_nnlo_as_0118_1000.yml: https://github.com/NNPDF/nnpdf/tree/master/n3fit/runcards/reproduce_nnpdf40/NNPDF40_nnlo_as_0118_1000.yml


Hessian conversion, compression, bundled sets and flavor number variations
--------------------------------------------------------------------------------
Among the released PDF sets of NNPDF4.0 are also some sets that are the result
of a transformation or combination of the fits produced using the
"fit runcards", these PDF grids are:

.. code-block::

    # Baseline NNDPF4.0 sets:
    NNPDF40_nnlo_as_0118
    NNPDF40_nnlo_as_0118_hessian

    # PDF sets with :math:`\alpha_s` variations:
    NNPDF40_nnlo_pdfas
    NNPDF40_nnlo_hessian_pdfas

    # PDF sets with flavor-number variations:
    NNPDF40_nlo_as_0118_nf_4
    NNPDF40_nlo_as_0118_nf_6
    NNPDF40_nnlo_as_0118_nf_4
    NNPDF40_nnlo_as_0118_nf_6
    NNPDF40_nlo_pch_as_0118_nf_3
    NNPDF40_nnlo_pch_as_0118_nf_3
    NNPDF40_nlo_as_0118_nf_4_pdfas
    NNPDF40_nnlo_as_0118_nf_4_pdfas

Section 10 of the NNPDF4.0 states how these can be obtained from the fits
produced using the method discussion under `Producing the NNPDF4.0 fits`_, but
we will again give some pointers here.

Both ``NNPDF40_nnlo_as_0118`` and ``NNPDF40_nnlo_as_0118_hessian`` are based on
a 1000 replica PDF set ``NNPDF40_nnlo_as_0118_1000``. Specifically,
``NNPDF40_nnlo_as_0118`` is the result of a compression of
``NNPDF40_nnlo_as_0118_1000`` using the
`pycompressor <https://github.com/N3PDF/pycompressor>`_ package, while
``NNPDF40_nnlo_as_0118_hessian`` can be created by running

.. code:: bash

    validphys NNPDF40_nnlo_as_0118_hessian.yml

For more information, see
:ref:`the tutorial on how to transform a Monte Carlo PDF set into a Hessian PDF set <mc2hessian>`.

The bundled PDF + :math:`\alpha_s` variation set ``NNPDF40_nnlo_pdfas`` can be
generated using the runcard |NNPDF40_nnlo_pdfas.yml|_, again, for more
information on how to bundle PDFs with :math:`\alpha_s` replicas, see
:ref:`the relevant tutorial <bundled-sets>`.


.. |NNPDF40_nnlo_pdfas.yml| replace:: ``NNPDF40_nnlo_as_0118_1000.yml``
.. _NNPDF40_nnlo_pdfas.yml: https://github.com/NNPDF/nnpdf/tree/master/n3fit/runcards/reproduce_nnpdf40/NNPDF40_nnlo_pdfas.yml


The PDF sets released as part of NNPDF4.0 also includes sets in which the 
maximum value of ``nf`` differs from the baseline value of ``nf=5``. To produce
these sets, the steps described in :ref:`howto nf variations`. In particular, 
the ID's of the required theories can be found in the 
`theory dabase <https://github.com/NNPDF/nnpdf/blob/master/nnpdfcpp/data/theory.db>`_,
where the id's of the theories used to do the flavor-number variations are 
218-227. 

.. warning::
    Please note that the Hessian and ``_pdfas`` conversion runcards won't work
    with the 4.0.3 version of the code (and thus the fixed enviroment) and
    instead require a version 4.0.4 or newer.

