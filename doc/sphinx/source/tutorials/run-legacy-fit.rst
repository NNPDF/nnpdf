.. _nnfit-usage:

How to run a legacy PDF fit (NNPDF3.1 style)
============================================

This tutorial explains how to run a PDF fit using the legacy code,
:code:`nnfit`. To find out how to run a PDF fit with the code that is currently
being developed, :code:`n3fit`, see the :ref:`n3fit-usage`.

0. Install the code, either using :ref:`conda<conda>` or by installing it from
:ref:`source<source>`. The former installation method should be considered the
default, with the latter being employed for situations in which you want to edit
the code yourself.

1. Create a runcard by, for example, using one of the configuration files
available `here <https://github.com/NNPDF/nnpdf/tree/master/nnpdfcpp/config/>`_.

2. Prepare the fit: use the command :code:`vp-setupfit --legacy <runcard>.yaml` to
generate a :code:`<runcard_folder>` in the current directory. This will have the
same name as the runcard and it will also contain a copy of the runcard itself.
The resources that are required to launch the fit, such as the theory ID and the
:math:`t_0` PDF set, will be downloaded automatically. Alternatively, such
resources can be obtained with the :ref:`vp-get<vp-get>` tool.

.. warning::

   The flat ``dataset_inputs`` structure for listing datasets in the runcard is
   not supported by ``nnfit``. Instead, datasets must be grouped explicitly by
   experiment, e.g.
   
.. code:: yaml

    experiments:
      - experiment: NMC
        datasets:
        - dataset: NMCPD
        - dataset: NMC
      - experiment: SLAC
        datasets:
        - dataset: SLACP
        - dataset: SLACD
   
3. To run the fit, use the :code:`nnfit` program. This takes a
:code:`<runcard_folder>` as input, as well a number that indexes the PDF replica
that will be produced. It is used as follows: :code:`nnfit <replica_number>
<runcard_folder>`. Therefore, to produce a fit with :math:`N` replicas, you will
need to run this command :math:`N` times with :code:`<replica_number>` running
from 1 to :math:`N`. Having to run this command :math:`N` times is ideal for
running a fit on a cluster.

4. Once each of the :code:`nnfit` commands has finished running, use
:code:`postfit <number_of_replicas> <runcard_folder>` to finalize the PDF set by
applying the post selection criteria. This will produce a set of
:code:`<number_of_replicas> + 1` replicas, since the mean of the fitted
replicas, usually dubbed :code:`replica_0`, will also be computed. Note that the
standard behaviour of :code:`postfit` can be modified by using various flags.
More information can be found at :ref:`postfit`.

5. You can then upload the fit to the
`server <https://data.nnpdf.science/fits/>`_ by using
:code:`vp-upload <runcard_folder>`. This will enable someone else to download
the fit using :ref:`vp-get<vp-get>`.

6. Finally, the fit can be analyzed with :code:`validphys`. You can create your
own custom analysis by referring to the :ref:`validphys guide<vp-index>`. To
produce a quick summary of the fit, you can compare it to a baseline fit by
using the :ref:`vp-comparefits <compare-fits>` command.
