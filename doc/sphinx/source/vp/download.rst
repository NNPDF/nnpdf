.. code:: {eval-rst}

   .. _download:

Downloading resources
=====================

``validphys`` is designed so that, by default, resources stored in known
remote locations are downloaded automatically and seamlessly used where
necessary. Available resources include PDF sets, completed fits,
theories, and results of past ``validphys`` runs that have been
`uploaded to the server <upload>`__. The ``vp-get`` tool, `described
below <#the-vp-get-tool>`__, can be used to download the same items
manually.

Automatic operation
-------------------

By default when some resource such as a PDF is required by ``validphys``
(or derived tools such as ``vp-setupfit``), the code will first look for
it in some local directory specified in the `profile
file <nnprofile>`__. If it is not found there, it will try to download
it from some remote repository (also specified in the profile).

For example a ``validphys`` runcard such as

.. code:: yaml

   pdf: NNPDF40_nnlo_as_01180
   fit: NNPDF40_nlo_as_01180

   theoryid: 208

   use_cuts: "fromfit"

   dataset_input:
       dataset: ATLASWZRAP36PB
       cfac: [EWK]

   actions_:
     - plot_fancy
     - plot_chi2dist

Will download if necessary the fit called ``NNPDF40_nlo_as_01180``, the
PDF set called ``NNPDF40_nnlo_as_01180`` and the theory with ID 208,
when validphys is executed with the default settings. In practice one
rarely has to worry about installing resources by hand when working with
NNPDF tools.

The behaviour of downloading automatically can be disabled by passing
the ``--no-net`` flag to supported tools. In that case, failure to find
a given resource locally will result in an error and exiting the
program. The ``--net`` flag makes the default behaviour explicit and has
no effect otherwise.

What can be downloaded
----------------------

The following resources are found automatically:

.. code:: {eval-rst}

   Fits
       Fits (specified by the ``fit`` key) can be downloaded if they have previously
       been uploaded with :ref:`vp-upload <upload>`. The corresponding PDF
       set will be installed as appropriate.

   PDF sets
       PDF sets (specified among others by the ``pdf`` key) are searched for in
       both NNPDF and LHAPDF repositories. If the PDF is not found and a fit with
       the same name exists, it will be downloaded and the corresponding PDF set
       will be installed and made available for usage.

   Theories
       Theories (specified by the ``theoryid`` key) are downloaded and
       uncompressed.

   ``validphys`` output files
       Files produced by ``validphys`` can be used as input to subsequent validphys
       analyses (for example χ² tables are used for αs fits). The user needs to
       have HTTP access to the repository, which is provided when installing using
       the :ref:`bootstrap script <conda>`. Output files are not specified by any
       top level config key, but instead actions can specify their own logic, for
       example for using an existing file instead of computing it.

.. code:: {eval-rst}

   .. _vp-get:

The ``vp-get`` tool
-------------------

The ``vp-get`` tool can be used to download resources manually, in the
same way ``validphys`` would do.

The basic syntax is

.. code:: bash

   vp-get <resource_type> <resource_name>

The available options for ``<resource type>`` can be seen with
``vp-get --list``. They correspond to the resources described
`above <#what-can-be-downloaded>`__.

.. code:: bash

   $ vp-get --list
   Available resource types:
    - fit
    - pdf
    - theoryID
    - vp_output_file

For example to download the fit ``NNPDF31_nlo_as_0118_1000`` we would
write

.. code:: bash

   $ vp-get fit NNPDF31_nlo_as_0118_1000

If the resource is already installed locally, the tool will display some
information on it and bail out:

.. code:: bash

   $ vp-get fit NNPDF31_nlo_as_0118_1000
   FitSpec(name='NNPDF31_nlo_as_0118_1000', path=PosixPath('/home/zah/anaconda3/envs/nnpdf-dev/share/NNPDF/results/NNPDF31_nlo_as_0118_1000'))

Downloading resources in code (``validphys.loader``)
----------------------------------------------------

.. code:: {eval-rst}

   The automatic download logic is implemented in the :py:mod:`validphys.loader`,
   specifically by the :py:class:`validphys.loader.RemoteLoader` and
   :py:class:`validphys.loader.FallbackLoader` classes.

   The logic is as follows: Given a resource type ``<foo>``, the normal
   :py:class:`validphys.loader.Loader` class would implement a ``check_<foo>`` method
   returning an object containing the appropriate metadata (such as file paths), or
   raise a ``LoaderError`` if the object cannot be found. The ``check_<foo>`` method
   of ``FallbackLoader`` (which is generated dynamically) will intercept the
   ``LoaderError`` and, if it happens, call the ``download_<foo>`` method of
   ``RemoteLoader``, if it exists. That method should cause the resource to be
   installed in such a way that the subsequent call of the ``Loader.check_<foo>``
   method succeeds. That is it should downoad the resource to the relevant search
   path, and uncompress it if needed.

In practice one can get a download aware loader by using a
``FallbackLoader`` instance, which will try to obtain all the required
resources from remote locations.

.. code:: python

   from validphys.loader import FallbackLoader as Loader

   l = Loader()
   #Will download theory 151 if needed.
   l.check_dataset('NMC', theoryid=151)

Conversely the ``Loader`` class will only search locally.

.. code:: python


   from validphys.loader import Loader

   l = Loader()

   l.check_dataset('NMC', theoryid=151)
   ---------------------------------------------------------------------------
   TheoryNotFound                            Traceback (most recent call last)
   <ipython-input-7-30e29a1539e8> in <module>
   ----> 1 l.check_dataset('NMC', theoryid=151)

   ~/nngit/nnpdf/validphys2/src/validphys/loader.py in check_dataset(self, name, rules, sysnum, theoryid, cfac, frac, cuts, use_fitcommondata, fit, weight)
       416 
       417         if not isinstance(theoryid, TheoryIDSpec):
   --> 418             theoryid = self.check_theoryID(theoryid)
       419 
       420         theoryno, _ = theoryid

   ~/nngit/nnpdf/validphys2/src/validphys/loader.py in check_theoryID(self, theoryID)
       288         if not theopath.exists():
       289             raise TheoryNotFound(("Could not find theory %s. "
   --> 290                   "Folder '%s' not found") % (theoryID, theopath) )
       291         return TheoryIDSpec(theoryID, theopath)
       292 

   TheoryNotFound: Could not find theory 151. Folder '/home/zah/anaconda3/share/NNPDF/data/theory_151' not found

Output files uploaded to the ``validphys`` can be retrieved specifying
their path (starting from the report ID). They will be either downloaded
(when using ``FallbackLoader``) or retrieved from the cache:

.. code:: python

   from validphys.loader import FallbackLoader as Loader
   l = Loader()
   l.check_vp_output_file('qTpvLZLwS924oAsmpMzhFw==/figures/f_ns0_fitunderlyinglaw_plot_closure_pdf_histograms_0.pdf')
   PosixPath('/home/zah/anaconda3/share/NNPDF/vp-cache/qTpvLZLwS924oAsmpMzhFw==/figures/f_ns0_fitunderlyinglaw_plot_closure_pdf_histograms_0.pdf')
