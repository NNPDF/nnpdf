.. _run-qed-fit:

==========================
How to run a QED fit
==========================
This tutorial describes how to run a QED fit using the LuxQED approach, as
described in `arXiv:1607.04266 <https://arxiv.org/abs/1607.04266>`_ and
`arXiv:1708.01256 <https://arxiv.org/abs/1708.01256>`_.

Setting up the runcard
----------------------

It is possible to perform a QED fit by adding a ``fiatlux`` block to the
runcard. The following code snippet shows an example of a QED fit configuration:

.. code-block:: yaml

    fiatlux:
      luxset: 251127-jcm-lh-qed-001
      additional_errors: true
      luxseed: 1234567890

The parameters contained in the ``fiatlux`` block are:

* ``luxset``
      The name of the PDF set used to generate the photon PDF with `FiatLux
      <https://github.com/scarrazza/fiatlux/>`_. The code will use as many
      photon replicas as the number of replicas contained in the ``luxset``. If
      the user tries to generate a replica with ID higher than the maximum ID of
      the ``luxset``, the code will start reusing photon replica from the first.
      Being the LuxQED approach an iterated procedure, and given that some
      replicas do not pass the ``postfit`` selection criteria, the user should
      make sure that the number of replicas in the ``luxset`` is high enough
      such that in the final iteration there will be a number of replicas higher
      than the final replicas desired.

* ``additional_errors``
      Boolean flag to switch on and off the additional errors of the LuxQED approach.

      .. note::

        The ``additional_errors`` flag should be switched on only in the very last
        iteration of the procedure.

* ``luxseed``
      The seed used to generate the additional errors of the LuxQED as in ``additional_errors``.

The code uses both the ``fiatlux`` block and the ``theoryid`` specified in the
runcard to identify the photon PDF set. As explained below, the code searches
for precomputed photon PDF sets using the pair of ``luxset`` and ``theoryid``
parameters, first locally and then on the NNPDF server (see
:ref:`photon-resources` for details). The photon PDF set will be named
``photon_theoryID_<theoryid>_fit_<luxset>``.

.. _photon-resources:

Running with Photon PDF sets
-----------------------------

The generation of a photon PDF set can add significant overhead to the fitting
procedure. Moreover, the same photon PDF set might be used in different fits. To
minimize the overhead due to the generation of photon PDF sets and avoid
redundant computations, the code looks for precomputed photon resources either
locally or on the NNPDF server. If a desired photon PDF set does not exist in
either of the two locations, it will be computed on the fly and stored locally.
The following sections describe how to use existing photon PDF sets or generate
new ones.

Using available Photon PDF sets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Before running a QED fit, it is strongly advised to prepare the fit using
``vp-setupfit``, as explained in :ref:`this tutorial <run-n3fit-fit>`. This will
ensure that all the required resources are available, including the photon PDF.
The desired photon PDF is specified by the ``luxset`` parameter in the
``fiatlux`` block and the ``theoryid``, as explained above. The code will first
look for the photon PDF set in the local directory specified in the
:ref:`profile file <nnprofile>`. If the set is not found locally, it will try to
download it from the NNPDF server. The following is an example of running
``vp-setupfit`` using the ``fiatlux`` block shown above:

.. code-block:: bash

    $ vp-setupfit qed_example_runcard.yml
    [INFO]: Could not find a resource (photonQED): Could not find Photon QED set /user/envs/nnpdf/share/NNPDF/photons_qed/photon_theoryID_702_fit_251127-jcm-lh-qed-001 in theory: 702. Attempting to download it.
    [INFO]: Downloading https://data.nnpdf.science/photons/photon_theoryID_702_fit_251127-jcm-lh-qed-001.tar.
    [==================================================] (100%)
    [INFO]: Extracting archive to /opt/homebrew/Caskroom/miniconda/base/envs/nnpdf/share/NNPDF/photons_qed
    [INFO]: Photon QED set found for 702 with luxset 251127-jcm-lh-qed-001.
    [WARNING]: Using q2min from runcard
    [WARNING]: Using w2min from runcard
    Using Keras backend
    [INFO]: All requirements processed and checked successfully. Executing actions.
    ...

This will download and extract the photon PDF set in the local
``photons_qed_path`` specified in the :ref:`profile file <nnprofile>`.

The ``vp-list`` utility tool can be used to list all the available photon PDF
sets locally and on the NNPDF server. To list the available photon PDF sets,
just run:

.. code-block:: bash

   $ vp-list photons
   [INFO]: The following photons are available locally:
    - theoryID_702_fit_251127-jcm-lh-qed-001
   [INFO]: The following photons are downloadable:
    - theoryID_702_fit_251127-jcm-lh-qed-002

In this example, there is one photon PDF set available locally, and one photon
resource available on the NNPDF server precomputed with theory ID 702 and
``251127-jcm-lh-qed-002`` as input PDF. The user can manually download a photon
PDF set using the ``vp-get`` tool as explained :ref:`here <download>`. For
example:

.. code-block:: bash

  $ vp-get photonQED 702 251127-jcm-lh-qed-002
  [INFO]: Could not find a resource (photonQED): Could not find Photon QED set /user/envs/nnpdf/share/NNPDF/photons_qed/photon_theoryID_702_fit_251127-jcm-lh-qed-002 in theory: 702. Attempting to download it.
  [INFO]: Downloading https://data.nnpdf.science/photons/photon_theoryID_702_fit_251127-jcm-lh-qed-002.tar.
  [==================================================] (100%)
  [INFO]: Extracting archive to /user/envs/nnpdf/share/NNPDF/photons_qed
  PosixPath('/user/envs/nnpdf/share/NNPDF/photons_qed/photon_theoryID_702_fit_251127-jcm-lh-qed-002')

As in the case of ``vp-setupfit``, this will download and extract the photon PDF
set in the local ``photons_qed_path`` specified in the :ref:`profile file <nnprofile>`.
Once the photon PDF set is available locally, the user can proceed to run the
fit with ``n3fit`` as usual.

.. warning::

   If ``vp-setupfit`` is not run before ``n3fit``, and the photon PDF set is not
   available locally, the code will **not** attempt to download it from the server,
   but will directly proceed to compute it on the fly. See the :ref:`next section <generating-photon-sets>` for
   more details.

.. _generating-photon-sets:
Generating new Photon PDF sets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If the desired photon PDF set is not available locally nor on the NNPDF server, the
code will generate the required photon PDF set replicas on the fly using `FiatLux <https://github.com/scarrazza/fiatlux/>`_.
This can be done either during the ``vp-setupfit`` step, which precomputes all photon
replicas before starting the fit, or during the fitting step with ``n3fit``, which
generates the photon replica as needed for each replica being fitted.

.. note::

  Generating photon PDF sets on the fly can add significant overhead to the
  fitting procedure. It is therefore strongly advised to precompute the photon
  PDF set using ``vp-setupfit`` before running the fit with ``n3fit``.

In either case, the generated photon PDF set will be stored locally in the
``photons_qed_path`` specified in the :ref:`profile file <nnprofile>`, so that
it can be reused in future fits. The folder containing the photon replicas will
be named as ``photon_theoryID_<theoryid>_fit_<luxset>`` where ``<theoryid>`` and
``<luxset>`` are the values specified in the runcard. The folder contains a ``npz``
file for each photon replica, named as ``replica_<replica_id>.npz``. Each replica
file contains the photon PDF grid computed with FiatLux at :math:`Q_{\rm init} = 100` GeV,
prior to the evolution through EKO.

.. important::

  Automatic upload to the NNPDF server through ``vp-upload`` is **not**
  supported at the moment. The user should manually create a ``tar`` archive
  file containing the photon replicas and upload it to server. Refer to the
  :ref:`profile file <nnprofile>` to find the remote URL where photon PDF sets are stored.


Using ``vp-setupfit`` (preferred)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  In order to trigger the computation of the photon PDF set during ``vp-setupfit``,
  the user needs to add the flag ``compute_in_setupfit: true`` in the ``fiatlux`` block
  discussed above. The following is an example of running ``vp-setupfit`` with
  this flag enabled:

  .. code-block:: bash

      $ vp-setupfit qed_example_runcard.yml
      [INFO]: Could not find a resource (photonQED): Could not find Photon QED set /user/envs/nnpdf/share/NNPDF/photons_qed/photon_theoryID_703_fit_251127-jcm-lh-qed-001 in theory: 703. Attempting to download it.
      [ERROR]: Resource not in the remote repository: Photon QED set for TheoryID 703 and luxset 251127-jcm-lh-qed-001 is not available in the remote server.
      [WARNING]: Photon QED set for theory 703 with luxset 251127-jcm-lh-qed-001 not found. It will be computed in vp-setupfit.
      [INFO]: Forcing photon computation with FiatLux during setupfit.
      [WARNING]: Using q2min from runcard
      [WARNING]: Using w2min from runcard
      Using Keras backend
      [INFO]: All requirements processed and checked successfully. Executing actions.
      ...

  In addition, the user can also specify the optional parameter ``eps_base`` to
  the ``fiatlux`` block to set the base precision of the FiatLux computation,
  which controls the precision on final integration of double integral. By
  default, it is set to ``1e-5``. If the ``compute_in_setupfit`` flag is set to
  ``true`` despite the photon PDF being available locally, the code will
  recompute and overwrite the existing photon PDF set

  .. tip::

    During ``vp-setupfit``, the code tries to distribute the computation of the
    photon replicas among the available CPU cores as specified in the
    ``maxcores`` key of the runcard. This can significantly speed up the
    computation of the photon PDF set. Make sure to set ``maxcores`` to a value
    compatible with the available hardware. For example, using ``maxcores: 64``
    will try to compute up to 64 photon replicas in parallel using 64 GB of RAM.

  Once the preparation step is completed, and the photon PDF set is generated and stored
  locally, the user can proceed to run the fit with ``n3fit`` as usual.

Using ``n3fit`` (discouraged)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  If the user prefers to compute the photon PDF set during the fitting step with
  ``n3fit``, no additional flag is needed in the runcard and ``vp-setupfit`` is
  not required beforehand (unless other resources are needed such as the
  :ref:`theory covariance matrix <vptheorycov-index>`). The code will check for
  the presence of the photon PDF set locally before starting the fit for each
  replica. If it is not found, it will proceed to compute the photon replica as
  needed for each replica being fitted. The following is an example of running
  ``n3fit`` where the photon PDF set is computed during the fitting step:

  .. code-block:: bash

      $ n3fit qed_example_runcard.yml 1
      [INFO]: Creating replica output folder in /user/runcards/qed_example_runcard/nnfit/replica_1
      [WARNING]: Using q2min from runcard
      [WARNING]: Using w2min from runcard
      Using Keras backend
      [WARNING]: No Photon QED set found for Theory 703 with luxset 251127-jcm-lh-qed-001. It will be computed using FiatLux. This may impact performance. It is recommended to precompute the photon set before running the fit. Refer to https://docs.nnpdf.science/tutorials/run-qed-fit.html for more details on precomputing photon PDF sets.
      [INFO]: All requirements processed and checked successfully. Executing actions.

  .. warning::

    Computing the photon PDF set during `n3fit` with multiple replicas or using GPUs
    has not been tested and may lead to unexpected behaviour. It is strongly advised to
    precompute the photon PDF set using ``vp-setupfit`` before running the fit with ``n3fit``.
