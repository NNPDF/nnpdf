.. _collect:

The collect function
====================

In the user interface we have the possibility to perform a computation
looping over a list of namespaces. In the code, we can define
providers that collect the results of such computations with the
`collect` function.

The signature is:

.. code:: python

	collect('provider_name', fuzzyspec)

This will expand the :ref:`fuzzyspec <fuzzyspecs>` relative to the current
:ref:`namespace <namespaces>` and
compute the function once for each frame.  Then it will put all the
results in a list (to be iterated in the same order as the fuzzyspec)
and set that list as the result of the provider. The provider in the
first argument is found following the standard `reportengine` rules.
It can be a function defined in a provider module, a configuration
input or a production rule, as well as another `collect` provider. As
a special case, one can pass functions directly
(defined with the `def` keyword). For example

.. code:: python

	possets_predictions = collect(positivity_predictions, ('posdatasets',))

Compared to a simple `for` loop, the collect function has the
advantages that the computations are appropriately reused and several
results could be computed simultaneously in the parallel mode.

We can use the output of `collect` as input to other providers. For
example:

.. code:: python

	def count_negative_points(possets_predictions):
	    return np.sum([(r.rawdata < 0).sum(axis=1) for r in
		    possets_predictions], axis=0)

`collect` can be used to appropriately group nested inputs. For
example, here is how to obtain a list of the experiments for each fit.

.. code:: python

	fits_experiments = collect('experiments', ('fits',))


Note that `collect` always returns a flat list with the provider
evaluated for each of the namespaces spanned by the fuzzyspec. For
example

.. code:: python

	fits_experiments_chi2_flat = collect(abs_chi2_data_experiment,
	    ('fits', 'fitcontext', 'experiments',))

results in a flat list
containing the result of `abs_chi2_data_experiment` resolved for each
experiment in each fit.  One can instead retain the structure by
chaining several `collect` providers. For instance, the code

.. code:: python

	experiments_chi2 = collect(abs_chi2_data_experiment, ('experiments',))
	fits_experiment_chi2_data = collect('experiments_chi2', ('fits', 'fitcontext'))

will result in `fits_experiment_chi2_data` producing one result for
each fit, where each of them is itself a list where each item result of
`abs_chi2_data_experiment` is evaluated for each experiment in a given
fit.

Standard iteration techniques can be used to process the results of
collect. For example, here is how we would print the χ² for each
experiment in each fit:

.. code:: python

	def print_fits_experiments_chi2(
		fits, fits_experiments, fits_experiment_chi2_data):
	    for fit, fit_experiments, experiments_data in zip(
		    fits, fits_experiments, fits_experiment_chi2_data):
		 print(f"Printing results for {fit}")
		 for experiment, chi2data in zip(fit_experiments, experiments_data):
		     print(f"χ² for {experiment} is ",
		        f"{chi2data.central_result}/{chi2data.ndata}")


A minimal runcard to use the action above is:

.. code:: yaml

	fits:
	  - NNPDF40_nlo_as_01180
	  - NNPDF40_nnlo_as_01180

	use_cuts: "fromfit"

	actions_:
	  - print_fits_experiments_chi2
