Getting started with validphys
==============================

To use ``validphys`` you must provide a YAML input runcard which includes

* The resources you need (PDFs, fits, etc.)
* The actions (functions) you would like to be carried out
* Additional flags and parameters for fine-tuning
* Metadata describing the author, title and keywords

To get an idea of the layout, :ref:`vpexamples` details the example runcards that can be found in
`this folder <https://github.com/NNPDF/nnpdf/tree/master/validphys2/examples>`_. The :ref:`tutorials`
section also takes you through how to make runcards for various tasks.

A simple example is:

.. code:: yaml

	pdf: NNPDF40_nlo_as_01180

	theoryid: 208

	use_cuts: "internal"

	dataset_input:
	    dataset: ATLASWZRAP36PB
	    cfac: [EWK]

	actions_:
	  - plot_fancy
	  - plot_chi2dist

We are specifying one PDF (by the `LHAPDF <https://lhapdf.hepforge.org/>`_ id),
one dataset and one
theory. Note that the dataset specification is identical to that of
the ``n3fit`` configuration files.

We are saying that we do not want to use any cuts on the data
(so we do not have to specify a fit containing the cut data, for example).

The special ``actions_`` key is used to declare the actions we want to
have executed. The syntax is the same as for the targets inside the
report (see :ref:`tut_report`).  We want a data-theory comparison (``plot_fancy``;
see :ref:`datthcomp`) and to
plot the distribution of the chi² for each replica (``plot_chi2dist``).

Once you have created a runcard (e.g. ``runcard.yaml``), simply run

.. code::

   validphys runcard.yaml

to set the ball rolling. For information on writing more complex runcards see
:ref:`here <complex_runcards>`.

Another useful command to be aware of is ``vp-comparefits - i``, which launches an interactive
session to compare two fits. See the tutorial :ref:`compare-fits` for more information.

For more tailored analysis, the API provides a high level interface to the code, allowing you to
extract objects and play around with them. See :ref:`vpapi`.

Finally, the ``validphys --help`` command can give you information on modules and specific actions, e.g.

.. code::

       $ validphys --help plots

will list all the actions defined in the plots module together with a brief description of each of them.
Asking for the help of one of the actions will list all the inputs required for this action. For example:

.. code::

   	$ validphys --help fits_chi2_table

   	fits_chi2_table

	Defined in: validphys.results

	Generates: table

	fits_chi2_table(fits_total_chi2_data, fits_datasets_chi2_table,
	  fits_groups_chi2_table, show_total: bool = False)

	Show the chi² of each and number of points of each dataset and
	experiment of each fit, where experiment is a group of datasets
	according to the ``experiment`` key in the PLOTTING info file, computed
	with the theory corresponding to the fit. Dataset that are not
	included in some fit appear as ``NaN``



	The following additionl arguments can be used to control the
	behaviour. They are set by default to sensible values:

	  show_total(bool) = False
	  per_point_data(bool) = True [Used by fits_groups_chi2_table]

We can see which keys have a special meaning in the configuration file with

.. code::

    $ validphys --help config

All other keys are interpreted literally (although they could be further processed by specific actions).
