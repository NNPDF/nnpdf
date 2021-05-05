.. _addvpplots:

How to add a new type of plot or table in ``validphys``
=======================================================

So you want to make a plot or table of PDFs, data, covariance matrices etc.?

Pre-existing actions
--------------------

**First check whether there is already an action which makes this type of plot.**

Some places to look:

-  `Here for plots/tables of data <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/dataplots.py>`_,
   and see :ref:`datthcomp` for guidance.
-  `Here for plots of PDFs <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/pdfplots.py>`_,
   and see the tutorial on plotting PDFs for guidance.
-  `Here for covariance matrix plots <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/theorycovariance/output.py>`_,
   and see :ref:`theory-covmat-examples` for guidance.
-  `Here for theory covariance validation plots <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/theorycovariance/tests.py>`_,
   and see :ref:`theory-covmat-examples` for guidance.
-  Use ``validphys --help``. This will give information on installed modules and also specific actions. For example,
   ``validphys --help plot_chi2_dist`` gives you information about the action ``plot_chi2_dist`` including input
   resources, output type and the place where it is defined.
     
Outline of approach
-------------------
**If there is no action for the plot you want, you can add your own!**

.. warning::
    Please follow the :ref:`rules`! In particular, this means you should create a new branch of the code
    to make your edits on, then create a pull request, with reviewers, for merging into the master branch.
    This helps to avoid bugged actions being propagated.
    
1. Locate the right place to put your plot/table. See the "places to look" above to find similar plots.
2. Identify the inputs you need for your plot/table. For example, if you are plotting some \\(\\chi^2\\)s you 
   may need the action ``chi2_data`` as an input. Try searching in the relevant modules in ``validphys``, using
   ``validphys --help`` or looking at the table of common inputs below.
3. See if you can make use of functions in `plotutils <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/plotutils.py>`_
   to help with making a plot.
4. Write an action which takes the inputs and returns a figure or table.
5. See :ref:`pytools` for general help on developing the python code. One approach is to make use of 
   the `IPython embed <https://ipython.readthedocs.io/en/stable/interactive/reference.html#embedding>`_
   function. This allows you to play around with variables from inside
   the action you are making. Another is to use an IPython environment such as a `Jupyter notebook <https://jupyter.org/>`_, 
   and obtain the necessary inputs using the :ref:`API <vpapi>`.
6. Use the ``@figure`` decorator for figures and ``@table`` decorator for tables. This means that
   the output will be saved in ``output/figures`` and ``output/tables`` respectively.

Common inputs
-------------

.. csv-table:: 
   :header: "Input name", "Description"
   :widths: 50, 50
   
   "``pdf`` or ``pdfs``", "PDF(s)"
   "``fit`` or ``fits``", "Fit(s)"
   "``processed_metadata_group``", "Name of the grouping for data, e.g. ``experiment`` or ``nnpdf31_process``"
   "``groups_index``", "Dataframe index for grouped data"
   "``groups_results``", "(Data, Theory) tuples for each group of data"
   "``groups_data``", "``DataSetSpec`` objects for each group"
   "``groups_central_values``", "Theory central values"
   "``groups_data_values``", "Data central values"
   "``groups_covmat``", "Experimental covariance matrix"
   "``theory_covmat_custom``", "Scale variation theory covariance matrix"
   "``total_chi2_data``", "``Chi2Data`` namedtuple with entries ``replica_result``, ``central_result`` and ``ndata`` for total \\(\\chi^2\\), **not including any correlations outwith experiments**"
   "``groups_chi2``", "``Chi2Data`` namedtuple for all groups"
   "``replica_data``", "PDF replicas"
   "``experiments_xq2map``", "Two ( \\(x\\), \\(Q^2\\) ) tuples, one for fitted data and one for cut data, for each of the ``dataset_inputs``"

Worked example for plotting \\(\\phi\\)
---------------------------------------
Say we want to make a plot which shows the breakdown of the \\(\\phi\\) values by group. This action actually already exists,
and it's called ``plot_phi``, but we will use it as an example. First we have a look in 
`dataplots <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/dataplots.py>`_ to see if there are similar 
actions. We can see that there is already an action called ``plot_groups_data_chi2`` which does something similar but for 
\\(\\chi^2\\)s instead. It looks like this:

.. code:: python

	@figure
	def plot_groups_data_chi2(groups_data, groups_chi2, processed_metadata_group):
	    """Plot the chi² of all groups of datasets with bars."""
	    exchi2 = []
	    xticks = []
	    for group, group_res in zip(groups_data, groups_chi2):
		exchi2.append(group_res.central_result/group_res.ndata)
		xticks.append(group.name)
	    fig, ax = plotutils.barplot(exchi2, collabels=xticks, datalabels=[r'$\chi^2$'])
	    ax.set_title(r"$\chi^2$ distribution by {}".format(processed_metadata_group))
	    return fig

To adapt this for \\(\\phi\\), we need to swap out the action ``groups_chi2`` for ``groups_data_phi`` instead. So the
start of the action will look like:

.. code:: python

	@figure
	def plot_phi(groups_data, groups_data_phi, processed_metadata_group):
	    """
	    We need to fill in this docstring later.
	    """
	    from IPython import embed
	    embed()

Note that an ``IPython embed`` has been included to help explore the various objects. This is for development purposes and we 
will remove it when making any commits. We can set up a test runcard that calls our action, like:

.. code:: yaml

	meta:
	    title: I didn't change the title
	    keywords: [Guilty]
	    author: Lazy Person

	pdf: {id: "NNPDF31_nnlo_as_0118", label: "3.1 NNLO"}

	theoryid: 53

	use_cuts : nocuts

	dataset_inputs:
	  - { dataset: NMC }
	  - { dataset: ATLASTTBARTOT, cfac: [QCD] }
	  - { dataset: CMSZDIFF12, cfac: [QCD,NRM], sys: 10 }  
	  - { dataset: LHCBWZMU7TEV, cfac: [NRM] }
	  - { dataset: LHCBWZMU8TEV, cfac: [NRM] }
	  - { dataset: ATLASWZRAP36PB}

	actions_:
	 - plot_phi


Running this runcard and using ``embed`` we can see that ``groups_data_phi`` is a list of tuples, one for each group. 
Each tuple is (``phi``,``num_points``). So we can write our action in a similar way to the ``plot_groups_data_chi2``
one:

.. code:: python

	@figure
	def plot_phi(groups_data, groups_data_phi, processed_metadata_group):
	    """plots phi for each group of data as a bar for a single
	    PDF input

	    See `phi_data` for information on how phi is calculated
	    """
	    phi = [exp_phi for (exp_phi, npoints) in groups_data_phi]
	    xticks = [group.name for group in groups_data]
	    fig, ax = plotutils.barplot(phi, collabels=xticks, datalabels=[r'$\phi$'])
	    ax.set_title(r"$\phi$ by {}".format(processed_metadata_group))
	    return fig
	    	  
Note that:

-  We used a list comprehension rather than a ``for`` loop.
-  We added a docstring explaining what the action does.
-  The action makes use of ``barplot`` from ``plotutils`` to do a bar chart.
-  The ``@figure`` decorator causes a figure to be generated and saved.	   
