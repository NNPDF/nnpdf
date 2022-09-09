.. _complex_runcards:

Writing `validphys` runcards
============================

In this section we go into some more detail on how to write `validphys`
runcards, in particular for more complex cases.

.. note::

	:ref:`vpexamples` details the example runcards that can be found in
	`this folder <https://github.com/NNPDF/nnpdf/tree/master/validphys2/examples>`_.
	The :ref:`tutorials` section also takes you through how to make runcards for
	various tasks.

We start with the following simple example:

.. code:: yaml

	pdf: NNPDF40_nnlo_as_01180

	theoryid: 208

	use_cuts: "internal"

	dataset_input:
	    dataset: ATLASWZRAP36PB
	    cfac: [EWK]

	actions_:
	  - plot_fancy
	  - plot_chi2dist

Multiple inputs and namespaces
------------------------------

Resources can be declared:

- At top level, like in the simple runcard above;
- Inside a mapping (with an arbitrary key);
- Inside an element of a list of mappings.

These mappings are called `namespaces`. For detailed information see :ref:`namespaces`.

.. important::
	When choosing your arbitrary key, good practice is to use a capital letter at the start.
	This helps to differentiate user-defined namespaces from internal objects.

1. Arbitrary namespaces

In this case we can modify the example as follows:

.. code:: yaml

	pdf: NNPDF40_nnlo_as_01180

	theoryid: 208

	fit: NNPDF40_nlo_as_01180

	With_cuts:
	  use_cuts: "fromfit"

	Without_cuts:
	  use_cuts: "nocuts"

	dataset_input:
	    dataset: ATLASWZRAP36PB
	    cfac: [EWK]

	actions_:
	  - With_cuts plot_fancy
	  - Without_cuts plot_chi2dist



Here `With_cuts` and `Without_cuts` are *arbitrary* strings that
specify *namespaces*.
We are asking for

- `plot_fancy` to be executed taking into account the cuts (note that we also need to
  specify the fit where they are read from)
- `plot_chi2dist` to be executed without the cuts.

Similar to
a programming language like C, the inner namespace has priority with
respect to the outer. For example, if we add a PDF specification to the
`with_cuts` namespace like this:


.. code:: yaml

	pdf: NNPDF40_nnlo_as_01180

	theoryid: 208

	fit: NNPDF40_nlo_as_01180

	With_cuts:
	  use_cuts: "fromfit"
	  pdf: NNPDF40_example_closure_test

	Without_cuts:
	  use_cuts: "nocuts"

	dataset_input:
	    dataset: ATLASWZRAP36PB
	    cfac: [EWK]

	actions_:
	  - With_cuts plot_fancy
	  - Without_cuts plot_chi2dist


The `plot_fancy` action will ignore the outer pdf
(NNPDF40\_nnlo\_as\_01180) and use the one defined in the innermost
namespace (NNPDF40_example_closure_test). Because we have not specified `plot_chi2dist` to
be executed within the `With_cuts` namespace, it will continue to use
NNPDF40\_nlo\_as\_01180.


2. Lists of namespaces

We can also have lists of mappings acting as namespaces. The action
will then be repeated inside each of the namespaces generating one
result for each. For example:

.. code:: yaml

	pdf: NNPDF40_nlo_as_01180

	theoryid: 208

	fit: NNPDF40_example_closure_test

	Specifications:
	- use_cuts: "fromfit"
	  pdf: NNPDF40_nnlo_as_01180

	- use_cuts: "nocuts"

	dataset_input:
	    dataset: ATLASWZRAP36PB
	    cfac: [EWK]

	actions_:
	  - Specifications plot_fancy

Now a different `plot_fancy` action will be executed for each of the
two mappings of the list "*Specifications*": one will use the NNLO PDF
and use the cuts from NNPDF40_example_closure_test, and the other will plot all points
in the dataset.

Some keys are appropriately interpreted either as lists of objects or
list or namespaces depending on the context. They are documented in
`validphys --help config`. For example, the `pdfs` key is entered as
a list of LHAPDF ids:

.. code:: yaml

	pdfs:
	  - NNPDF40_nlo_as_01180
	  - NNPDF40_nnlo_as_01180


Because the `plot_fancy` action takes a list of pdfs as input,
something like this:

.. code:: yaml

	pdfs:
	  - NNPDF40_nlo_as_01180
	  - NNPDF40_nnlo_as_01180

	theoryid: 208

	use_cuts: "nocuts"

	dataset_input:
	    dataset: ATLASWZRAP36PB
	    cfac: [EWK]

	actions_:
	  - plot_fancy


will produce plots where the two PDFs appear together. However,
we can also produce individual plots for each PDF, by simply
specifying that we want to loop over `pdfs`:

.. code:: yaml

	pdfs:
	  - NNPDF40_nlo_as_01180
	  - NNPDF40_nnlo_as_01180

	theoryid: 208

	use_cuts: "nocuts"

	dataset_input:
	    dataset: ATLASWZRAP36PB
	    cfac: [EWK]

	actions_:
	  - pdfs plot_fancy


In this case the value of the `pdfs` key is seen as equivalent to:

.. code:: yaml

	pdfs:
	  - {pdf: NNPDF40_nlo_as_01180}
	  - {pdf: NNPDF40_nnlo_as_01180}


However, the special treatment allows us to simplify both the input
file and the programmatic interface of the functions.

Nesting namespaces
------------------

Namespace specifications like those described above can be arbitrarily
nested. Values will be searched from the inner to the outer namespace. When
the namespace specifications represent lists of mappings, all possible
combinations will be produced.

Consider the example:

.. code:: yaml

	pdfs:
	    - NNPDF40_nlo_as_01180
	    - NNPDF40_nnlo_as_01180
	    - NNPDF40_nnlo_as_01180_hessian

	fit: NNPDF40_nlo_as_01180

	theoryids:
	    - 208
	    - 162

	With_cuts:
	    use_cuts : "nocuts"

	dataset_inputs:
	    - { dataset: LHCBWZMU7TEV, cfac: [NRM] }
	    - { dataset: LHCBWZMU8TEV, cfac: [NRM] }
	    - { dataset: ATLASWZRAP36PB }

	actions_:
	  - With_cuts::theoryids::pdfs::dataset_inputs plot_fancy

This will first enter the "*With_cuts*" namespace (thus setting
``use_cuts = "nocuts"`` for the action), and then loop over all the
theories, pdfs and datasets.

The order over which the looping is done is significant:

1. The outer specifications must set all the variables required for the inner
   ones to be fully resolved (so `With_cuts` must go before `dataset_inputs`).

2. The caching mechanism works by grouping together the namespace
   specifications from the beginning. For example, suppose we were to
   add another action to the example above:

.. code:: yaml

    - with_cuts:
        theoryids:
          pdfs:
            dataset_inputs:
                - plot_chi2dist

both of these require the same convolutions to be computed. `Validphys` will
realize this as long as both actions are iterated in the same way.
However, permuting `pdfs` and `theoryids` would result in the
convolutions computed twice, since the code cannot prove that they
would be identical.

In summary:
 - Always loop from more general to more specific.
 - Always loop in the same way.

Action arguments
----------------

Action arguments are syntactic sugar for specifying arguments visible
to a single action. They are subject to being verified by the action-defined
checks. For example, in the PDF plotting example above:

.. code:: yaml

	pdfs:
	    - NNPDF40_nlo_as_01180
	    - NNPDF40_nnlo_as_01180
	    - NNPDF40_nnlo_as_01180_hessian

	First:
	    Q: 1
	    flavours: [up, down, gluon]

	Second:
	    Q: 100
	    xgrid: linear

	actions_:
	  - First::plot_pdfreplicas (normalize_to=NNPDF40_nlo_as_01180)
	  - First plot_pdfs
	  - Second plot_pdfreplicas


The `normalize_to` key only affects the `plot_pdfreplicas` action.
Note that defining it inside the `first` mapping would have had the
same effect in this case.


The `from_` special key
-----------------------

The `from_` special key specifies that the value of a resource is to be taken from
a container. This is useful for working with fits (but not limited to
that). For example:

.. code:: yaml

	fit: NNPDF40_nlo_as_01180

	use_cuts: "nocuts"

	description:
	    from_: fit

	theory:
	    from_: fit

	theoryid:
	    from_: theory

	Q: 10

	template: report.md

	normalize:
	    normalize_to: 1

	datanorm:
	    normalize_to: data

	pdfs:
	    - from_: fit
	    - NNPDF40_nnlo_as_01180

	data_inputs:
	    from_: fit

	actions_:
	   - report(out_filename=index.md)


Here the `from_` key is used multiple times:

 - To obtain the description string from the report input card.
 - To obtain the theory mapping from the fit input card.
 - To obtain the theoryid key from the theory mapping.
 - To obtain a single PDF produced in the fit (as an element of the
   list/namespaces of pdfs). Note that the keyword is also allowed
   inside nested elements.
 - To obtain a set of all the experiments of the fit.

The `from_` key respects lazy processing, and therefore something like
this will do what you expect:

.. code::  yaml

	fits:
	    - NNPDF40_nlo_as_01180
	    - NNPDF40_nnlo_lowprecision

	use_cuts: "nocuts"

	theory:
	    from_: fit

	theoryid:
	    from_: theory

	Q: 10

	description:
	    from_: fit

	template: report.md

	normalize:
	    normalize_to: 1

	datanorm:
	    normalize_to: data

	pdfs:
	    - from_: fit
	    - NNPDF40_nnlo_as_01180_hessian

	dataset_inputs:
	    from_: fit

	actions_:
	  - fits report

This will work exactly as the example above, except that a new action
(with its corresponding different set of resources) will be generated
for each of the two fits.

For fits, there is a shortcut to set `dataset_inputs`, `pdf` and
`theoryid` to the values obtained from the fit. This can be done with
the `fitcontext` rule. The above example can be simplified like this:

.. code:: yaml

	fits:
	    - NNPDF40_nlo_as_01180
	    - NNPDF40_nnlo_lowprecision

	use_cuts: "nocuts"

	Q: 10

	description:
	    from_: fit

	template: report.md

	normalize:
	    normalize_to: 1

	datanorm:
	    normalize_to: data

	pdfs:
	    - from_: fit
	    - NNPDF40_nnlo_as_01180_hessian

	actions_:
	  - fits::fitcontext report

Note that one still needs to set manually other keys like `description` and `pdfs`.

``from_: Null``
--------------

As a special case, `from_: Null` will retrieve the variable from the
current namespace. This comes handy to transform lists of items into
other items. Consider for example:

.. code:: yaml

	Base:
	    fit: NNPDF40_nnlo_as_01180_1000

	Pairs:
	    fits:
		- from_: Base
		- from_: null

	fits:
	    - NNPDF40_nnlo_as_01180_NNPDF31
	    - NNPDF40_nnlo_as_01180_collider_only
	    - NNPDF40_nnlo_as_01180_DIS_only
	    - NNPDF40_nnlo_as_01180_nojets
	    - NNPDF40_nnlo_as_01180_noLHCbb
	    - NNPDF40_nnlo_as_01180_noLHC
	    - NNPDF40_nnlo_as_01180_notop
	    - NNPDF40_nnlo_as_01180_noZpT
	    - NNPDF40_nnlo_as_01180_nophoton
	    - NNPDF40_nnlo_as_01180_ATLASW8TeV
	    - NNPDF40_nnlo_as_01180_noATLASCMSDY
	    - NNPDF40_nnlo_as_01180_EMC

	use_cuts: "fromfit"

	printopts:
	    print_common: False

	description:
	    from_: fit

	meta:
	    author: Zahari Kassabov
	    keywords: [nn40final, gallery]

	template_text: |
	    % Non-default datasets

	    The datasets are compared to the default `{@Base fit@}` fit.

	    {@with fits::fitcontext@}
	    {@fit@}
	    ======

	    {@description@}

	    {@with Pairs@}

	    {@printopts print_dataset_differences  @}
	    {@print_different_cuts@}

	    {@endwith@}
	    {@endwith@}

	actions_:
	  - report(main=True, mathjax=True)

- At the beginning, we are printing the name of the fit contained in
  `Base`.
- Then we are iterating over each of the `fits` (that we
  defined explicitly in the config), and using `fitcontext` to set some
  variables inside the `with` block.
- In the inner block `{@with Pairs@}`, we are making use of the definition
  of `Pairs` to set the `fits` variable to contain two fits: the one defined in `Base` and the
  one that changes with each iteration.
- Because the actions `print_dataset_differences` and `print_different_cuts` are inside that
  `with` block, the value of the variable `fits` they see is precisely
  this pair, which supersedes our original definition, inside that
  block.

The `namespaces_` special key
-----------------------------

The `namespaces_` key can be used to form a list of namespaces in
a similar way as with the `{@with@}` block in the report. A key difference
is that the `namespaces_`
block allows the list to be names, and in this way it can interact
with providers expecting a complex input structure. The namespace
elements are separated by `::` and have the same meaning as in the
report.  Consider the following example:

.. code:: yaml

	dataspec_input:
	  - fitdeclarations:
	       - NNPDF40_nlo_as_01180
	       - NNPDF40_nnlo_as_01180
	    fits_computed_psedorreplicas_chi2_output: new-alldata/fits_matched_pseudorreplicas_chi2_table.csv
	    fits_chi2_paramfits_output: new-alldata/central_global.csv
	    badspecs:
	       - badcurves: discard
		 speclabel: "Global, discard"
	       - badcurves: allminimum
		 speclabel: "Global, allminimum"

	  - fitdeclarations:
	       - NNPDF31_nnlo_as_0117_uncorr_collider
	       - NNPDF31_nnlo_as_0118_uncorr_collider
	    fits_computed_psedorreplicas_chi2_output: new-alldata/collider.csv
	    fits_chi2_paramfits_output: new-alldata/collider_central.csv
	    badspecs:
	       - badcurves: discard
		 speclabel: "Collider, discard"
	       - badcurves: allminimum
		 speclabel: "Collider, allminimum"

	dataspecs:
	    namespaces_: "dataspec_input::badspecs
		          ::fits_as_from_fitdeclarations::fits_name_from_fitdeclarations
		          ::use_fits_computed_psedorreplicas_chi2_output::use_fits_chi2_paramfits_output"

	meta:
	   author: Zahari Kassabov
	   title: Summary of the allminimum and discard for global and collider only fits
	   keywords: [as]

	template_text: |

	    We compare the results of the determinations with `allminimum`
	    and `discard` on the global and collider only fits.

	    # Table

	    {@dataspecs_as_value_error_table@}

	    # Plot

	    {@plot_dataspecs_as_value_error@}

	actions_:
	  - report(main=True)



Here we are generating a list of namespaces called `dataspecs` which
the actions `dataspecs_as_value_error_table` and
`plot_dataspecs_as_value_error` expect as an input, starting from the
product of each of the two elements in the `dataspec_input` list and
its corresponding `badspecs` inner namespace, so that we have four
namespaces in total, labelled "Global, discard", "Global, allminimum",
"Collider, discard" and "Collider, allminimum". We are further
applying production rules  to extract the
information we need from the fit names and input files, producing the
corresponding values inside the correct `dataspecs` entry.

The whole list namespace is then passed as input to the actions (which
are implemented using :ref:`the collect function <collect>`).

This advanced functionality allows us to generate almost arbitrary
inputs in a declarative way and using very few primitives, at the cost
of a bit of learning curvature.

Currently the `namespaces_` functionality is restricted to generating
namespaces that are used at top level.

Plotting labels
---------------

Several resources (PDFs, theories, fits) support a short form where
one specifies the ID required to recover the resource (e.g. LHAPDF ID,
theory ID and fit folder respectively) and also form where a plotting
layer is specified together with the ID. For example:

.. code:: yaml

	pdfs:
	    - id:  NNPDF40_nlo_as_01180
	      label: NLO

	    - id: NNPDF40_nnlo_as_01180
	      label: NNLO

	    - id: NNPDF40_nnlo_as_01180_hessian
	      label: Hessian NNLO


In all plots the label will be used everywhere the PDF name needs to
be displayed (like in legends and axes).

The plotting labels for datasets are read from the `dataset_label` key
in the plotting files.

See :ref:`pdfplots` for examples.
