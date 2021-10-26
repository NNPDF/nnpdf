 .. _design:

The design of ``validphys 2``
=============================

Introduction
------------

The objective of this page is to explain the high level considerations
underpinning the design of ``validphys 2``, and of ``reportengine``
:cite:p:`zahari_kassabov_2019_2571601`, the code it is based on. It should be
useful for anyone aiming to understand the general philosophy and design goals
of the project or take in in qualitatively different directions.  More concrete
information on how to use the code can be found in the sections under Using
validphys.

Some specific issues of scientific code
---------------------------------------

There are some characteristics of typical scientific programs that are
rather specific to them. A program aiming to optimize the workflow of a
scientific collaboration should be designed with those in mind. Some of
these characteristics make work easier compared to other types of
programs like text editors or games. In general the runtime execution
can be made simpler because user interactivity is not required. On the
other hand we need easy access to all parts of the code, since we never
know which parts might be subject to changes. That’s the whole point of
research!

We now proceed to list these characteristics in more detail.

Lack of interactivity
~~~~~~~~~~~~~~~~~~~~~

Typically programs such as report generators, fitting codes, or
simulations are entirely determined by the initial user input together
with perhaps some random source and environment events like which
replicas failed. We typically don’t make network connections in the
middle of fits, don’t wait for the user to be pushing buttons. In many
instances it is possible to know about the amount of resources such as
memory that a scientific program is going to require shortly after
running it.

Changing structure
~~~~~~~~~~~~~~~~~~

It happens often that the hardcoded parameter you are optimizing becomes
the next object of research. Scientific code should make as little
assumptions as possible on how it is supposed to run, while managing to
get by and produce some results. This somewhat upsets the balance of
trade-offs compared to other kinds of programs. We are willing to trade
some increased complexity for the ability to look under the hood.

Complexity
~~~~~~~~~~

Non trivial scientific programs can rarely be assembled exclusively by
putting together external components without adding anything of our own.
That typically consists on unspecified, or at best, underspecified
algorithms that are difficult to test in isolation.

Design considerations
---------------------

Look before you leap
~~~~~~~~~~~~~~~~~~~~

A scientific program usually requires a large set of complex input
parameters and can take a very long time to produce results. Results
tend to be frequently unexpected even if everything is working
correctly. It is essential to check that the inputs are correct and make
sense. Thus good error reporting capabilities are an important design
consideration.

Also, the inputs of a given function should be checked as early as
possible, which is not necessarily at the point of the program where the
function is to be called. For example, something like this:

.. code:: python

   def check_parameter_is_correct(parameter):
       ...

   def plot(complex_calculation, parameter):
       check_parameter_is_correct(parameter)
       #make plot
       ...

has the disadvantage that in order to get to the check, we presumably
need to compute ``complex_calculation`` first, and that could be a waste
of time if it turns out that ``parameter`` is not correct for some
reason and we can’t do the plot anyway. Instead we should arrange to
call ``check_parameter_is_correct(parameter)`` as early as possible
(outside the plotting function) and show the user an informative error
message.

   All inputs should be checked as early as possible. If the checks
   pass, the program should not fail.

Of course it is impossible to guarantee that this is always going to
happen. For example it is possible that we check that a folder exists
and we have permissions to write to it, but by the time we actually need
to write, the folder has been deleted. However in such cases the state
of the program can be considered broken and it’s OK to make it fail
completely.

There is an API for early checks in ``reportengine``. We would write
something like:

.. code:: python

   @make_argcheck
   def check_parameter_is_correct(parameter):
       ...

   @check_parameter_is_correct
   def plot(complex_calculation, parameter):
       #make plot
       ...

The checking function will now be called as soon as the program realizes
that the plot function will be required eventually (at “compile time”).
The user would be shown immediately an error message explaining why the
parameter is wrong.

The fancy term for this style of coding is `Contract
Programming <https://en.wikipedia.org/wiki/Design_by_contract>`__.

Declarative input
~~~~~~~~~~~~~~~~~

It is convenient to be able to specify the *what* the program should do
without any regard of knowledge of *how* that is achieved by the
underlying implementation. The primary input of validphys are
`YAML <https://en.wikipedia.org/wiki/YAML>`__ run cards. A very simple
one looks like this:

.. code:: yaml

    pdfs:
        - NNPDF40_nlo_as_01180
        - NNPDF40_nnlo_as_01180
        - NNPDF40_nnlo_as_01180_hessian

   norm:
       normalize_to: NNPDF40_nlo_as_01180

   first:
       Q: 1
       flavours: [up, down, gluon]

   second:
       Q: 100
       scale: linear

   actions_:
       - first::norm plot_pdfreplicas
       - first plot_pdfs
       - second plot_pdfreplicas

This has a number of advantages:

Correct by definition
^^^^^^^^^^^^^^^^^^^^^

A declarative input specifies what you want. It is up to the underlying
code to try to provide it (or fail with an informative message).

Clear meaning
^^^^^^^^^^^^^

It is easy for a human to verify that the input is indeed what it was
intended. Even without any explanation it should be easy enough to guess
what the runcard above does.

Implementation independent
^^^^^^^^^^^^^^^^^^^^^^^^^^

The input is very loosely coupled with the underlying implementation,
and therefore it is likely to remain valid even after big changes in the
code are made. For example, in the runcard above, we didn’t have to
concern ourselves with how LHAPDF grids are loaded, and how the values
of the PDFs are reused to produce the different plots. Therefore the
underlying mechanism could change easily without breaking the runcard.

Usage as a programmatic API
~~~~~~~~~~~~~~~~~~~~~~~~~~~

While the goal of ``reportengine`` is to allow simple and easily
repeatable bach actions, sometimes it is far simpler to get the work
done with a raw (Python) script, or it is needed to explore the outcomes
using something like an Jupyter notebook. It would be good to be able to
use all the tools that already exist in ``validphys`` for that, without
needing to reinvent the wheel or to alter functions so that for example
they don’t write data to some preconfigured path.

Therefore The various computing and plotting tools should work well when
included in a normal script that doesn’t use the ``reportengine`` graph
compiler.

This is implemented by making sure that as much as possible all the
``validphys`` functions are *pure*. That is, the output is a
deterministic function of the inputs, and the function has no side
effects (e.g. no global state of the program is altered, nothing is
written to disk). There are some exceptions to this though. For example
the function that produces a reweighted PDF set needs to write the
result to disk. The paths and side effects for other more common results
like figures are managed by ``reportengine``. For example, the
``@figure`` decorator applied to a function that returns a Python
(``matplotlib``) figure will make sure that the figure is saved in the
output path, with a nice filename, while having no effect at all outside
the ``reportengine`` loop. The same goes for the check functions
described above.

Easy to loop
~~~~~~~~~~~~

Very frequently there is the need to compare. It is easy enough to write
simple scripts that loop over the required configurations, but that
cannot scale well when requirements change rapidly (and is also easy to
make trivial mistakes). Therefore reportengine allows configurations to
easily loop over different sets of inputs. For example the following
runcard:

.. code:: yaml

   pdfs:
       - id:  NNPDF40_nlo_as_01180
         label: NLO

       - id: NNPDF40_nnlo_as_01180
         label: NNLO


   theoryids:
       - 208
       - 200
   use_cuts : nocuts

   experiments:
     - experiment: LHCb
       datasets:
         - { dataset: LHCBWZMU7TEV, cfac: [NRM] }
         - { dataset: LHCBWZMU8TEV, cfac: [NRM] }

     - experiment: ATLAS
       datasets:
         - { dataset: ATLASWZRAP36PB}

   actions_:
    - theoryids::pdfs::experiments::experiment plot_fancy

Will produce a separate plot for each combination of the two theories
(200 and 208), the two PDFs at the top, and each dataset in the two
experiments (so 28 plots in total). This syntax is discussed in more
detail in the [Usage] section.

It should be trivial to repeat an action for different sets of inputs.
