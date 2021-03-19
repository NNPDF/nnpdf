Defining custom pipelines
-------------------------

Here we discuss what needs to go from user entered strings in the YAML
file plots and reports.

The basic code flow is as follows:

 1. The `actions_` key is parsed to obtain a list of requirements with
	their associated [fuzzyspec](#fuzzyspecs).

 2. Each requirement spans other requirements. These can be:
    - Providers: Other functions with requirements on their own.
	- User input processed by the [Configuration], which is
	  immediately tested for correctness.
	- Production rules, also derived from the configuration.

 3. Once the requirements are satisfied for a given provider, the
	[checks](#checking-providers) of the provider are executed.

 4. If all the checks pass, all the runtime requirements are executed
	in such an order that the dependencies are resolved.

### Configuration

A configuration class derived from `reportengine.ConfigParser` is used to parse the
user input. In validphys, it is defined in `validphys.Config`.

The parsing in reportengine is *context dependent*. Because we want to
specify resources as much as possible before computing anything (at
"*compile time*"), we need to have some information about other
resources (e.g. theory folders) in order to do any meaningful
processing.

The `Config` class takes the user input and the dependencies and:

 - Returns a valid resource if the user input is valid.

 - Raises a `ConfigError` if the user input is invalid.

To parse a given user entered key (e.g. `posdataset`), simply define
a `parse_posdataset` function. The first argument (i.e. second after
`self`) will be the raw value in the configuration file. Any other
arguments correspond to dependencies that are already resolved at the
point where they are passed to the function (`reportengine` takes care
of that).

For example, we might have:
```python
def parse_posdataset(self, posset:dict, * ,theoryid):
    ...
```


The type specification (`:dict` above) makes sure that the user input
is of that type before it is seen by the function (which avoids
a bunch of repetitive error checking). A positivity dataset requires
a theory ID in order to be meaningfully processed (i.e. to find the
folder where the fktables are) and therefore the theoryid will be
looked for and processed first.

We need to document what the
resource we are providing does. The docstring will be seen in
`validphys --help config`:
```python
def parse_posdataset(self, posset:dict, * ,theoryid):
    """An observable used as positivity constrain in the fit.
    It is a mapping containing 'dataset' and 'poslambda'."""
    ...
```

#### Production rules

Apart from `parse_` functions, which take an explicit user input from
the corresponding key (and optionally a set of dependencies), there
are the `produce_` functions, which take only the dependencies. Other
than not taking the user input, the `produce_` functions work in
a very similar way to the `parse_` functions: They are resolved at
*"compile time"*, before any procider function is executed, and  they
should raise a `ConfigError` if they fail.

In general production rules should be preferred to parse functions
that bundle together various dependencies (e.g. data, cuts and
theory), because by having more granular elements, we can iterate over
them in different ways: For examples we might want to generate
a separate report page for each of the positivity datasets, where they
are compared for multiple theories. We could break the parse function
above into:

```python
def parse_posdataset_input(self, posset:dict):
    ...

def produce_posdataset(posdataset_input, *, theoryid):
   ...
```

Now the user has to enter a key called "posdataset_input", from which
some Python object will be obtained as the return value of
`parse_posdataset_inout`. Then, `produce_posdataset` is used to an
object representing the positivity set and the corresponding FKTables
in a given theory is obtained from the output of
`parse_posdataser_input` and a theory ID.

#### Automatically parsing lists

It is possible to easily process list of elements once the parsing for
a single element has been defined. Simply add an `eleement_of`
decorator to the parsing function defined in the Config class:
```python
@element_of('posdatasets')
def parse_posdataset(self, posset:dict, * ,theoryid):
```

Now `posdatasets` is parsed as a list of positivity datasets, which
can be passed together to a provider, or iterated over, (for example
with a `with` tag in the report, see [Report template specification]).

Note that you can also put together results from evaluating providers
using [the collect function], which can be used to map computations
over the lists described here.



#### Validphys loaders

In `validphys`, we use a `Loader` class to load resources from various
folders. It is good to have a common interface, since it is used to
list the available resources of a given type or even download
a missing resource. The functions of type `check_<resource>` should
take the information processed by the Config class anf verify that
a given resources is correct. If so they should return a "Resouce
specification" (something typically containing metadata information
such as paths, and a `load()` method to get the C++ object from
`libnnpdf`). We also define a `get` method that returns the C++ object
directly (although I am not sure it's very useful anymore).

In the case of the positivity set, this is entirely given in terms of
existing check functions:

```python
def check_posset(self, theiryID, setname, postlambda):
    cd = self.check_commondata(setname, 0)
    fk = self.check_fktable(theiryID, setname, [])
    th =  self.check_theoryID(theiryID)
    return PositivitySetSpec(cd, fk, postlambda, th)

def get_posset(self, theiryID, setname, postlambda):
    return self.check_posset(theiryID, setname, postlambda).load()
```

A more complicated example should raise the appropriate loader
errors (see the other examples in the class).

The `PositivytySetSpec` could be defined roughly like:
```python
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
```
Here `PositivitySet` is the `libnnpdf` object. It is generally better
to pass around the spec objects because they are lighter and have more
information (e.g. the theory in the above example).

With this, our parser method could look like this:
```python

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
```

The first part makes sure that the user input is of the expected form
(a mapping with a string and a number). The `ConfigError` has support
for suggesting that something could be mistyped. The syntax is
`ConfigError(message, bad_key, available_keys)`. For example, if the
user enters "poslanda" instead of "postlambda", the error message
would suggest the correct key.

Note that all possible error paths must end by raising
a `ConfigError`.



### Computing PDF dependent quantities

Now that we can receive positivity sets as input, let's do something
with them. The SWIG wrappers allow us to call the C++ methods of
`libnnpdf` from Python. These things go in the `validphys.results`
module. We can start by defining a class to produce and hold the
results:
```python
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
```

`pdf.stats_class` allows to interpret the results of the convolution
as a function of the PDF error type (e.g. to use the different
formulas for the uncertainty of Hessian and Monte Carlo sets). In that
way it allows to abstract away the different error types. One
constructs an object inheriting from `validphys.core.Stats` that is
appropriate for a given error type by calling `pdf.stats_class(data)`
where data is an array where the entries along the first dimension are
the results from each member computed from `libnnpdf` (and the other
dimensions are arbitrary). `Stats` has methods that appropriately
collapse along the first axis. For example `central_value` computes
the mean along the first axis for Monte Carlo PDFs and yields the
first member for Hesssian PDFs.

And then define a simple provider function:
```python
def positivity_predictions(pdf, positivityset):
     return PositivityResult.from_convolution(pdf, positivityset)
```

### The collect function

In the user interface we have the possibility to perform a computation
looping over a list of namespaces. In the code, we can define
providers that collect the results of such computations with the
`collect` function.

The signature is:
```python
collect('provider_name', fuzzyspec)
```

This will expand the `fuzzyspec` relative to the current namespace and
compute the function once for each frame.  Then it will put all the
results in a list (to be iterated in the same order as the fuzzyspec)
and set that list as the result of the provider. The provider in the
first argument is found following the standard `reportengine` rules.
It can be a function defined in a provider module, a configuration
input or a production rule, as well as another `collect` provider. As
a special case, one can pass directly functions
(defined with the `def` keyword).  For example
```python
possets_predictionsa = collect(positivity_predictions, ('posdatasets',))
```

Compared to a simple `for` loop, the collect function has the
advantages that the computations are appropriately reused and several
results could be computed simultaneously in the parallel mode.

We can use the output of `collect` as input to other providers. For
example:
```python
def count_negative_points(possets_predictions):
    return np.sum([(r.rawdata < 0).sum(axis=1) for r in
	    possets_predictions], axis=0)
```

`collect` can be used to appropriately group nested inputs. For
example here is how to obtain a list of the experiments for each fit.
```python
fits_experiments = collect('experiments', ('fits',))
```

Note that `collect` always returns a flat list with the provider
evaluated for each of the namespaces spanned by the fuzzyspec. For
example
```python
fits_experiments_chi2_flat = collect(abs_chi2_data_experiment,
    ('fits', 'fitcontext', 'experiments',))
```
results in a flat list
containing the result of `abs_chi2_data_experiment` resolved for each
experiment in each fit.  One can instead retain the structure by
chaining several `collect` providers. For instance, the code
```python
experiments_chi2 = collect(abs_chi2_data_experiment, ('experiments',))
fits_experiment_chi2_data = collect('experiments_chi2', ('fits', 'fitcontext'))
```
will result in `fits_experiment_chi2_data` producing one result for
each fit, where each of them is itself list where each item result of
`abs_chi2_data_experiment` evaluated for each experiment in a given
fit.

Standard iteration techniques can be used to process the results of
collect. For example here is how we would print the χ² for each
experiment in each fit:

```python
def print_fits_experiments_chi2(
        fits, fits_experiments, fits_experiment_chi2_data):
    for fit, fit_experiments, experiments_data in zip(
            fits, fits_experiments, fits_experiment_chi2_data):
         print(f"Printing results for {fit}")
         for experiment, chi2data in zip(fit_experiments, experiments_data):
             print(f"χ² for {experiment} is ",
                f"{chi2data.central_result}/{chi2data.ndata}")
```

A minimal runcard to use the action above is:

```yaml
fits:
  - NNPDF31_nlo_as_0118
  - NNPDF31_nnlo_as_0118

use_cuts: "fromfit"

actions_:
  - print_fits_experiments_chi2
```


### Checking providers

Providers can checks that verify that all the required preconditions
are met. Checks are executed at the time at which the call node is
just created and all its required dependencies are either in the
namespace or scheduled to be produced. Checking functions take the
current state of the namespace, as well as an unspecified set of other
parameters (because I haven't decided on the interface yet!).
Therefore check functions should accept `**kwargs` arguments. Checks
are decorated with the `reportengine.checks.make_argcheck` function.
If checks don't pass, they must raise
a `reportengine.checks.CheckError` exception.

For example, given a reweighting function, we may want to check that
the current PDF (the value that will be passed to the function) has
a Monte Carlo error type, we might define a check like:
```python
@make_check
def check_pdf_is_montecarlo(ns, **kwargs):
    pdf = ns['pdf']
    etype = pdf.ErrorType
    if etype != 'replicas':
        raise CheckError("Error type of PDF %s must be 'replicas' and not %s"
                          % (pdf, etype))

```

Checks can be used (abused) to modify the namespace before the action
function sees it. This can be used for some advanced context dependent
argument default setting (for example setting default file names based
on the nsspec).

The check is associated to the provider function by simply applying it
as a decorator:

```python
@check_pdf_is_montecarlo
def chi2_data_for_reweighting_experiments(pdf, args):
    ...
```

A slightly higher level interface to checks is implemented by the
`make_argcheck` decorator. Instead of receiving a namespace and other
unspecified arguments, like the functions decorated with `make_check`,
it simply takes the arguments we want to test. The function can return
a dictionary that will be used to update the namespace (but that is
not required, it can also not return anything).

For example the `check_pdf_is_montecarlo` above could be more easily
implemented like:
```python
@make_argcheck
def check_pdf_is_montecarlo(pdf):
    etype = pdf.ErrorType
    if etype != 'replicas':
        raise CheckError("Error type of PDF %s must be 'replicas' and not %s"
                          % (pdf, etype))
```
`make_argcheck` should be preferred, since it is more explicit, and
could be extended with more functionality later on. However it is
newer and not very used currently in the code.

Checks have no effect outside of reportengine (unless you call them
explicitly).

Ideally, the checks should be sufficient to guarantee that the
actions will not fail at runtime.

### Producing figures

In order to produce figures, just decorate your functions returning
`matplotlib` `Figure` objects  with the `reportengine.figure.figure`
function, e.g.:

```python
@figure
def plot_p_alpha(p_alpha_study):
   fig, ax = plt.subplots()
   #Plot something
   ...
   return fig
```

This will take care of the following:

 - Saving the figures with a nice, unique name to the output folder,
   in the formats specified by the user.

 - Closing the figures to save memory.

 - Making sure figures are properly displayed in reports.

There is also the `figuregen` decorator for providers that are
implemented as generators that yield several figures (see e.g. the
implementation of `plot_fancy`). Apart from just the figure, yield
a tuple (prefix, figure) where the prefix will be used in the
filename.

### Producing tables

These work similarly to [figures](#producing-figures), as described
above. Instead use the `@table` and `@tablegen` decorators.

Tables will be saved in the CSV formats.

### Customizing how things look in the report

By default, the `str()` method will be applied to objects that appear
in the report. If you want a custom behaviour, declare a declare
a custom `as_markdown` property for your objects. It should return
a string in Pandoc Markdown describing your object. Raw HTML is
also allowed (although that decreases the compatibility, e.g. if we
decide to output LaTeX instead of HTML in the future).
