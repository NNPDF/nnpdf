How to reproduce an NNPDF4.0 fit
================================================================================

Here we describe how to reproduce the NNPDF4.0 PDF sets that are publicly
available through LHAPDF as discussed in section 10 of the NNPDF4.0 paper
:cite:p:`nnpdf40`.


The files needed to reproduce the results are available in the folder
`n3fit/runcards/reproduce_nnpdf40
<https://github.com/NNPDF/nnpdf/tree/master/n3fit/runcards/reproduce_nnpdf40>`_
of the `project repository <https://github.com/NNPDF/nnpdf>`_ on Github. The
`.yml` files referenced on this page can be found in that folder.

To start, we need to set up a :ref:`conda environment <conda>` that is the same
as the one used to produce the NNPDF4.0 PDF sets, meaning it is an a conda
environment with the version of ``nnpdf`` and the versions of all its
dependencies are equal to those originally used to produce the NNPDF4.0 PDF
sets. Such an environment can be generated and activated by running

.. code:: bash

    conda env create --file nnpdf40_env.yaml
    conda activate nnpdf40

Next, using this environment, we can produce the fits in the
:ref:`usual way <run-n3fit-fit>`. The names of the runcards in the
``reproduce_nnpdf40/`` folder are the same as those of the those of the
corresponding (public) LHAPDF grids of the same name. For example, baseline fit
with :math:`\alpha_s(m_Z)=0.118` and a  variable-flavor-number scheme with up to
five active flavors can be generated with ``NNPDF40_nnlo_as_0118_1000.yml``.
