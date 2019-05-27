Report Engine
=============

The main components are:

 - A Report object. 
 - An interface to construct a report from a given configuration:
	 `ResourceProcessor`.
 - A convention for the layout of the code that actually produces the
	 report, based on `ResourceProcessor` plus some helper tools to enforce
	 that convention (`repottools`?).

The main interface should be coded in Python 3.5.

Here a graphical representation:
![PDF](https://github.com/NNPDF/reportengine/tree/designdocs/docs/reportengine_layout.pdf)

Report
------

The Report object should be as simple as possible and easy to extend.

The basic interface of the report object will consist on methods to:

- Extract the required resources and parameters (see below) from the
	user specification.

- Fill the construct the report with the resources provided by the
	engine.
Our implementation will use `jinja2` templates to compile to Markdown.
It will also be responsible for handling more complex objects, such as
`matplotlib` figures and tables (for example to export them to the
correct path).

The report class should only understand a set of basic objects and
know how to compile them to Markdown. This will call the appropriate
macros in `jinja2` (for example, place a caption below a figure).

Default helpers and styles exist to produce the final content (i.e.
PDF or html files out of Markdown) but are largely a concern of the
clients.

Together with the template, data file, which is processed by
a `ResourceProcessor` is used to supply the information necessary to
build the report.

ResourceProcessor
-----------------

This class takes a configuration file which contains *input resources*
and figures out the processing that needs to be done to obtain *output
resources* (which are specified in the report layout in this case). 

While it is made with reports in mind, it is way more general.

In addition to tools for running the configuration, there are tools to
validate the inputs and analyze the requirements.

While in the most simple case, the outputs are directly functions of
the inputs, there is also the possibility of intermediate steps, so
that the required operations can be represented as a *Direct Acyclic
Graph (DAG)*.

The class also takes a Python package which contains the functions
needed to actually carry out the computations (as well as doing
specific validations). This package will contain three kinds of
functions:

- Parsers: Parse the user input in the configuration file to produce
	an input resource.

- Checks: Perform domain specific checks on the configuration before
	it is executed.

- Providers: Compute a specific resource (taking other resources
	and parameters as input).

Each of these is described in more detail below.

###Example

Imagine we have one *input resource*, *fit*. And we need to produce
two *output resources*: *arclength plots* and *arclength affinty*.

Both of these require a rather expensive computation called
*arclength*.

`ResourceProcessor` would be called with these parameters, as well as
a module `validphys` containing the following *provider* functions:

```python

#provides the *resource* arclength
def arclength(fit):
    ...

def arclength_affinity(arclength):
    ...

def plot_arclength(arclength):
    ...

```

The first job of `ResourceProcessor` is to figure out the appropriate
order in which these functions need to be called. 


###Documentation

The docstrings of the provider functions (plus some additions) are
automatically made available as documentation in the various
interfaces. For example:

```python

#provides the *resource* arclength
def arclength(fit):
    """Computes the arclength for all replicas in ``fit``. It uses
	``nndiff`` to calculate the analytical expression of the
	derivative of the neural network.
    ...

def plot_arclength(arclength):
    """Plot distribution of arclengths."""
    ...
```

Would cause the docs to be automatically available to various parts in
the code and possibly the report itself (tooltips of figures?).


###Namespaces and parameters

It is possible to call *final* providers with different *parameters* in
different parts of the report.
For example, imagine a function:

```python
def plot_compare_pdfs(pdf, base_pdf):
    ...
```

It might be useful to call that function with `base_pdf` pointing at
the previous fit or at the closure test prior in different parts of
the report. As long as this resource is *final* (i.e. if the compare
PDF figures are not required by any other resource in the report).

All resources are resolved within *namespaces* which can be specified
by the user (in the report layout). They are used to produce resources
with completely different input, which would correspond to different
execution graphs, in general.

There is also the default namespace, which will be used if no
namespace is specified.

An independent DAG will be constructed for each
namespace and all the functions necessary to construct them will be
executed (even if the inputs are the same). Any caching behaviour is
responsibility of the client.

###Interactive interfaces

Using a `DAGResources` should be trivial to determine what needs to be
recomputed when some *input resource* is changed by the user. This
sets the stage for building a more interactive representation such as
a web interface.

###Checks

As much as possible failures due to incorrect user input must be
checked **before** any computation takes place. `ResourceProcessor` should
implement the basic checks (i.e. that the graph can actually be
constructed). Additionally the Python 3.5 [`typing
module`](https://docs.python.org/3/library/typing.html) could be used
to check for the types of the resources.

There is also support for domain-specific checks implemented as
a `check` decorator. It takes a function as a parameter which in turn
is called with the decorated function, the namespace and the instance
of `ResourceProcessor`. The decorated function  For example:

```python
def pdfs_installed(resource, namespace, resolver):
   """Check if the relevant pdfs are installed in LHAPDF"""
   ...

@check(pdfs_installed)
def plot_compare_pdfs(pdf:PDF, base_pdf:PDF) -> reportengine.Figure:
    ...
```

The fact that the arguments are in fact PDFs would be checked by
`ResourceProcessor` (which will know the return types of all producers),
while the function `pdfs_installed` would be called before actually
building the report.

The checks are called in the same order as the functions would.

###Input

The input resources are set by the user with a YAML file.  Keys ending
with an underscore *_* have a special meaning and are not allowed for
clients. One such key is `namespaces_`, which is used to declare
namespaces (see above).  A basic configuration file would be:

```yaml
fit: /path/to/fit
base_pdf: NNPDF30_as_118

namespaces_:
    vsprevious:
	    base_pdf: prevfitpdf
	vsclosure:
	    base_pdf: MMHT
```

This would create 3 namespaces (the two explicitly defined and the
default one), each of which contain two input
resources,  *fit* and *base_pdf*. Fit is equal for all of them, and of
inherited from the global namespace, while *base_pdf* is different for
each of them.

An application specific parser would be defined to process the
resource. It is implemented by the functions `parse_<resource>`
defined in the client package (maybe with the option of specifying
which module). They take as input the value only:

```python
from valifphys.core import Fit, PDF

def parse_fit(self, path) -> Fit:
	return Fit(path)

def parse_base_pdf(self, lhapdfname) -> PDF:
	return PDF(lhapdfname)

```

If no such function is found, the resources will have the value given
in the YAML file. The type would  be checked as described above.

Any error during the parsing would be intercepted and recast as a nice
error message.

###SMPDF correspondence

Many of these ideas are directly taken from
[SMPDF](https://github.com/scarrazza/smpdf). In particular the
[`actions`](https://github.com/scarrazza/smpdf/blob/master/src/smpdflib/actions.py)
module is a primitive implementation of `ResourceProcessor`, though much of
the work is done manually.

A rough correspondence in terminology would be:

action -> provider

actiongroup -> namespace

Configuration -> ResourceProcessor

Eventually that part of SMPDF would be reworked to use this framework.

Script conventions
------------------

The client functions (*resource providers*) will take a set of
resources as input, and produce one new resource as output. The name
of the new resource will be the same as that of the function. They
will not have access to other information such as the `ResourceProcessor`
instance that is running them.

If they specify the return type using function signatures, it will be
checked before executing the graph. However this is not mandatory.

The function defined by the user should have no side effects. In
particular they should not modify their input resources in any way
(so we executing the DAG in any valid order would produce the same
result) and the output should be a deterministic function of the
output. This largely confines the possible errors in the client
applications to the functions where they have originated.

Also, if this is fulfilled, things like parallel execution of the DAG,
or fault tolerance (producing a valid report even if a function fails)
becomes feasible.

In particular the functions should not write files (such as images) to
disk,but instead should rely on centrally provided functionality
(which would take care of saving to the correct paths, save in the
requested format and so on).

The output of the final resources should be objects from
`reportengine`, or objects that can be casted to them. For example
a `reportenfine.Figure` object would contain one or more `matplotlib`
figures, but also the corresponding `caption` attributes.  Saving the
actual files to disk will be responsibility of `reportengine`.

Of course the *no side effects* rule has exceptions. One example is
caching the outputs, as indexed by the inputs (as implemented by
`functools.lru_cache` for example). However the it should *look like*
no side effects exist to the external world.

In addition to resource providers, the client can also implement
**checks** (see above). A difference is that checks do have access to
the namespace and the `ResourceProcessor` instance where they belong. This
is necessary to implement less trivial checks. For instance SMPDF
checks that the PDF requested for plotting exist in LHAPDF. However if
a new PDF set is to be generated in some previous step (and in another
"namespace"), it is also considered valid to enter the name of that
new PDF.

It is not clear what is the best strategy for aggregating multiple
objects of the same kind (for example compare the arclength
distribution of two different PDF sets). The current approach is to
make the resource providers explicitly aware and always operate over
lists of items. However other strategies might be considered, such as
a more general design of the graph-based execution and checking model.
