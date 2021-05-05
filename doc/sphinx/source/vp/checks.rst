Checking providers
==================

Providers can have checks that verify that all the required preconditions
are met. Checks are executed at the time at which the call node has
just been created and all of its required dependencies are either in the
namespace or scheduled to be produced. Checking functions take the
current state of the namespace, as well as an unspecified set of other
parameters (the interface is as yet undecided!).
Therefore check functions should accept `**kwargs` arguments. Checks
are decorated with the `reportengine.checks.make_argcheck` function.
If checks do not pass, they must raise
a `reportengine.checks.CheckError` exception.

For example, given a reweighting function, we may want to check that
the current PDF (the value that will be passed to the function) has
a Monte Carlo error type. We might define a check like:

.. code:: python

	@make_check
	def check_pdf_is_montecarlo(ns, **kwargs):
	    pdf = ns['pdf']
	    etype = pdf.ErrorType
	    if etype != 'replicas':
		raise CheckError("Error type of PDF %s must be 'replicas' and not %s"
		                  % (pdf, etype))


Checks can be used (abused) to modify the namespace before the action
function sees it. This can be used for some advanced context dependent
argument default setting (for example setting default file names based
on the nsspec).

The check is associated to the provider function by simply applying it
as a decorator:

.. code:: python

	@check_pdf_is_montecarlo
	def chi2_data_for_reweighting_experiments(pdf, args):
	    ...


A slightly higher level interface to checks is implemented by the
`make_argcheck` decorator. Instead of receiving a namespace and other
unspecified arguments, like the functions decorated with `make_check`,
it simply takes the arguments we want to test. The function can return
a dictionary that will be used to update the namespace (but that is
not required, it can also not return anything).

For example, the `check_pdf_is_montecarlo` above could be more easily
implemented like:

.. code:: python

	@make_argcheck
	def check_pdf_is_montecarlo(pdf):
	    etype = pdf.ErrorType
	    if etype != 'replicas':
		raise CheckError("Error type of PDF %s must be 'replicas' and not %s"
		                  % (pdf, etype))

`make_argcheck` should be preferred, since it is more explicit, and
could be extended with more functionality later on. However, it is
newer and not currently used very much in the code.

Checks have no effect outside of reportengine (unless you call them
explicitly).

Ideally, the checks should be sufficient to guarantee that the
actions will not fail at runtime.
