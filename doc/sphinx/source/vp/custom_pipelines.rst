Defining custom pipelines
=========================

Here we discuss what information from user-entered strings needs to go in the YAML
file to plots and reports.

The basic code flow is as follows:

 1. The `actions_` key is parsed to obtain a list of requirements with
    their associated :ref:`fuzzyspec <fuzzyspecs>`

 2. Each requirement spans other requirements. These can be:

    - Providers: Other functions with requirements of their own.
    - User input processed by the configuration, which is immediately tested for correctness.
    - Production rules, also derived from the configuration.

 3. Once the requirements are satisfied for a given provider, the
    checks of the provider are executed.

 4. If all the checks pass, all the runtime requirements are executed
    in such an order that the dependencies are resolved.

Configuration
-------------

A configuration class derived from `reportengine.ConfigParser` is used to parse the
user input. In validphys, it is defined in
`validphys.config <https://github.com/NNPDF/nnpdf/blob/master/validphys2/src/validphys/config.py>`_.

The parsing in reportengine is *context dependent*. Because we want to
specify resources as much as possible before computing anything (at
"*compile time*"), we need to have some information about other
resources (e.g. theory folders) in order to do any meaningful
processing.

The `config` class takes the user input and the dependencies and:
 - Returns a valid resource if the user input is valid.
 - Raises a `ConfigError` if the user input is invalid.

To parse a given user-entered key (e.g. `posdataset`), simply define
a `parse_posdataset` function. The first argument (i.e. second after
`self`) will be the raw value in the configuration file. Any other
arguments correspond to dependencies that are already resolved at the
point where they are passed to the function (`reportengine` takes care
of that).

For example, we might have:

.. code:: python

	def parse_posdataset(self, posset:dict, * ,theoryid):
	    ...


The type specification (`dict` above) makes sure that the user input
is of that type before it is seen by the function (which avoids
a bunch of repetitive error checking). A positivity dataset requires
a theory ID in order to be meaningfully processed (i.e. to find the
folder where the :ref:`FK tables <fktables>` are) and therefore the `theoryid` will be
looked for and processed first.

We need to document what the
resource we are providing does. The docstring will be seen in
``validphys --help config``:

.. code:: python

	def parse_posdataset(self, posset:dict, * ,theoryid):
	    """An observable used as positivity constrain in the fit.
	    It is a mapping containing 'dataset' and 'poslambda'."""
	    ...


Production rules
-----------------

Apart from `parse_` functions, which take an explicit user input from
the corresponding key (and optionally a set of dependencies), there
are the `produce_` functions, which take only the dependencies. Other
than not taking the user input, the `produce_` functions work in
a very similar way to the `parse_` functions: they are resolved at
*"compile time"*, before any provider function is executed, and they
should raise a `ConfigError` if they fail.

In general, production rules should be preferred to parse functions
that bundle together various dependencies (e.g. data, cuts and
theory), because by having more granular elements, we can iterate over
them in different ways: for example, we might want to generate
a separate report page for each of the positivity datasets, where they
are compared for multiple theories. We could break the parse function
above into:

.. code:: python

	def parse_posdataset_input(self, posset:dict):
	    ...

	def produce_posdataset(posdataset_input, *, theoryid):
	   ...


Now the user has to enter a key called "posdataset_input", from which
some Python object will be obtained as the return value of
`parse_posdataset_input`. Then `produce_posdataset` is used for an
object representing the positivity set and the corresponding FK tables
in a given theory is obtained from the output of
`parse_posdataser_input` and a theory ID.

Automatically parsing lists
---------------------------

It is possible to easily process list of elements once the parsing for
a single element has been defined. Simply add an `element_of`
decorator to the parsing function defined in the Config class:

.. code:: python

	@element_of('posdatasets')
	def parse_posdataset(self, posset:dict, * ,theoryid):


Now `posdatasets` is parsed as a list of positivity datasets, which
can be passed together to a provider, or iterated over (for example
with a `with` tag in the report, see :ref:`reports`).

Note that you can also put together results from evaluating providers
using the collect function, which can be used to map computations
over the lists described here.

Validphys loaders
-----------------

In `validphys`, we use a `Loader` class to load resources from various
folders. It is good to have a common interface, since it is used to
list the available resources of a given type or even download
a missing resource. The functions of type `check_<resource>` should
take the information processed by the Config class and verify that
a given resource is correct. If so, they should return a "Resource
specification" (something typically containing metadata information
such as paths, and a `load()` method to get the C++ object from
`libnnpdf`). We also define a `get` method that returns the C++ object
directly.

In the case of the positivity set, this is entirely given in terms of
existing check functions:

.. code:: python

	def check_posset(self, theoryID, setname, postlambda):
	    cd = self.check_commondata(setname, 0)
	    fk = self.check_fktable(theiryID, setname, [])
	    th =  self.check_theoryID(theiryID)
	    return PositivitySetSpec(cd, fk, postlambda, th)

	def get_posset(self, theoryID, setname, postlambda):
	    return self.check_posset(theiryID, setname, postlambda).load()


A more complicated example should raise the appropriate loader
errors (see the other examples in the class).

The `PositivitySetSpec` could be defined roughly like:

.. code:: python

	 class PositivitySetSpec():
	     def __init__(self, commondataspec, fkspec, poslambda, thspec):
		 self.commondataspec = commondataspec
		 self.fkspec = fkspec
		 self.poslambda = poslambda
		 self.thspec = thspec

	     @property
	     def name(self):
		 return self.commondataspec.name

	     def __str__(self):
		 return self.name

	     @functools.lru_cache()
	     def load(self):
		 cd = self.commondataspec.load()
		 fk = self.fkspec.load()
		 return PositivitySet(cd, fk, self.poslambda)

Here `PositivitySet` is the `libnnpdf` object. It is generally better
to pass around the spec objects because they are lighter and have more
information (e.g. the theory in the above example).

With this, our parser method could look like this:

.. code:: python

	def parse_posdataset(self, posset:dict, * ,theoryid):
	    """An observable used as positivity constrain in the fit.
	    It is a mapping containing 'dataset' and 'poslambda'."""
	    bad_msg = ("posset must be a mapping with a name ('dataset') and "
		       "a float multiplier(poslambda)")

	    theoryno, theopath = theoryid
	    try:
		name = posset['dataset']
		poslambda = float(posset['poslambda'])
	    except KeyError as e:
		raise ConfigError(bad_msg, e.args[0], posset.keys()) from e
	    except ValueError as e:
		raise ConfigError(bad_msg) from e

	    try:
		return self.loader.check_posset(theoryno, name, poslambda)
	    except FileNotFoundError as e:
		raise ConfigError(e) from e


The first part makes sure that the user input is of the expected form
(a mapping with a string and a number). The `ConfigError` has support
for suggesting that something could be mistyped. The syntax is
`ConfigError(message, bad_key, available_keys)`. For example, if the
user enters "poslanda" instead of "postlambda", the error message
would suggest the correct key.

Note that all possible error paths must end by raising
a `ConfigError`.

Computing PDF-dependent quantities
----------------------------------

Now that we can receive positivity sets as input, let's do something
with them. The SWIG wrappers allow us to call the C++ methods of
`libnnpdf` from Python. These things go in the `validphys.results`
module. We can start by defining a class to produce and hold the
results:

.. code:: python

	class PositivityResult(StatsResult):
	    @classmethod
	    def from_convolution(cls, pdf, posset):
		loaded_pdf = pdf.load()
		loaded_pos = posset.load()
		data = loaded_pos.GetPredictions(loaded_pdf)
		stats = pdf.stats_class(data.T)
		return cls(stats)

	    @property
	    def rawdata(self):
		return self.stats.data


`pdf.stats_class` allows for the interpretation of the results of the convolution
as a function of the PDF error type (e.g. to use the different
formulas for the uncertainty of Hessian and Monte Carlo sets). In that
way it allows to abstract away the different error types. One
constructs an object inheriting from `validphys.core.Stats` that is
appropriate for a given error type by calling `pdf.stats_class(data)`,
where data is an array where the entries along the first dimension are
the results from each member computed from `libnnpdf` (and the other
dimensions are arbitrary). `Stats` has methods that appropriately
collapse along the first axis. For example, `central_value` computes
the mean along the first axis for Monte Carlo PDFs and yields the
first member for Hesssian PDFs.

And then define a simple provider function:

.. code:: python

	def positivity_predictions(pdf, positivityset):
	     return PositivityResult.from_convolution(pdf, positivityset)
