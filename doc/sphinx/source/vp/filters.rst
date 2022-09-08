.. code:: {eval-rst}

   .. _filters:

Filtering data
==============

Introduction
------------

In PDF fits, not all the data provided by the experimental
collaborations are useful. For example, we may wish to discard certain
datapoints for which we know small-x resummation or electroweak
corrections are important. These effects are problematic since we know
them to be important, but we cannot account for them.

In this light, we produce cuts of the data, by filtering data points
which we know are free of the above and other problems.

In ``validphys 2``, the cuts are handled by the ``validphys.filters``
alongside filter definitions and defaults found within
``validphys.cuts``. .

Cuts as declarative filters
---------------------------

Due to the nature of data cuts, it is important to be transparent about
which cuts are being applied to which dataset and/or process. Moreover,
it is useful for the rules defining the data cut to be readable such
that a non-developmental user can read and understand the nature of the
rule by making these rules functions of kinematic variables such as
``p_T`` or ``Q2``.

In much the same vein, it is useful for any default values used in the
rules to be readily accessible. For example, suppose there is a minimum
value for the square transferred momenta in the DIS process ``q2min``,
that is used widely by many different rules. It is important for this
variable to be in an obvious and easily accessed location.

Defaults
--------

There are certain values which are commonly used by many rules. For
example, the value ``q2min`` usually takes the value ``3.49`` or
``w2min`` is usually set to ``12.5``.

It is thus useful to define these default values somewhere. These values
can be found within ``validphys.cuts`` inside the ``defaults.yaml``
file. One can overwrite these values and this is discussed later.

Filters
-------

In ``validphys 2`` the default filter rules used can be found in the
``validphys.cuts`` module within the ``filter.yaml`` file. This file is
read by ``validphys`` and is interpreted as a ``list`` of
``dictionaries``.

By default, these filters can have several entries:

1. ``dataset``: The dataset this rule applies
2. ``process_type``: The process type this rule applies to
3. ``rule``: The ``Python`` code defining the rule for this filter
4. ``reason``: (optional) The reason this rule was needed
5. ``local_variables``: (optional) Any additional, non-standard local
   variables the user wishes to add for this rule only.

.. code:: {eval-rst}

   .. note::
     At least one of :code:`dataset` or :code:`process_type` is required.
     Additionally, a :code:`rule` entry is always required.

The ``rule`` entry in the rule definition is ``eval``\ uated as
``Python`` code. If the rule does not apply to this particular datapoint
(say the dataset names don’t match) then we return ``None`` indicating
this rule had nothing to do with this particular datapoint. In this
case, we move on to the next rule. However, if the process type or
dataset defined in the rule match that of the datapoint, we evaluate the
rule. If the rule evaluates to ``False`` we discard the point, if
instead it returns ``True`` we move on to the next rule. If by the time
all the rules have been evaluated and we have yet to return ``False``,
then the datapoint passes and it is kept.

In addition, the user can add any theory parameter they wish. For
example, one could add ``PTO: NNLO`` which means to evaluate the rule
only if the theory is NNLO. These are discussed further `here <#PTO>`__.
One can see a full list of possible theory parameters using:
``vp-checktheory <theory id>``

.. code:: {eval-rst}


   .. important::
       The :code:`rule` entry should be interpreted as a :code:`str` type within :code:`Python`. As such
       a rule such as :code:`rule: True` is not valid since this is read in as a boolean,
       however, :code:`rule: "True"` is perfectly valid notation. Moreover, the string
       itself should be valid :code:`Python` code.

By default the user can use the following non-builtin mathematical
functions in their rules: ``sqrt``, ``log`` or ``fabs`` (floating point
absolute value). In addition, one can use any ``numpy`` function using
``np.<function>`` in their rule definition. For example:

.. code:: yaml

   rule: "np.exp(x) > 0.1"

The kinematic variables that can be used within the rule depends on the
process type. A full list of available parameters can be found by
running:

.. code:: ipython

   In [1]: from NNPDF import CommonData                                               

   In [2]: print(dict(CommonData.kinLabel))

The user may additionally define their own variables by adding the
``local_variables`` field to their rule. For example, I can use ``w2``
in my rule, so long as I define what I mean by ``w2``:

.. code:: yaml

     local_variables:
       w2: Q2 * (1 - x) / x

.. code:: {eval-rst}

   .. danger::
     Defining :code:`local_variables` is non-commutative. The order of definition is important.
     If a local variable depends on other local variables, then the user must ensure all other
     dependencies have already been defined.

The following would raise an error

.. code:: yaml

     local_variables:
       w: sqrt(w2)
       w2: Q2 * (1 - x) / x

The following would not

.. code:: yaml

     local_variables:
       w2: Q2 * (1 - x) / x
       w: sqrt(w2)

.. code:: {eval-rst}

   .. note::
     :code:`local_variables` have a local scope. They apply to only the rule within which
     they are defined.

Theory parameters and perturbative orders
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are particular situations in which we only want to evaluate a rule
if the theory input for the PDF matches certain conditions. For example,
it may be the case we only keep the datapoint provided the theory
includes intrinsic charm or is evaluated at NNLO.

Suppose for example I wish the rule to only be evaluated if the theory
includes intrinsic charm. We note in the ``theory.get_description()``,
the relevant entry is ``'IC': 1`` (we use here theory 53 for
demonstration purposes). Thus if I want my rule to be applied only if
the theory has intrinsic charm, I simply add to my rule:

.. code:: yaml

     IC: True

Similarly I can condition on flavour number scheme. I again check
``theory.get_description()`` and note that the relevant ``key`` is
``'FNS'``. Thus to only evaluate my rule if the FNS is ``FONLL-C``,
simply add:

.. code:: yaml

     FNS: FONLL-C

Similarly, one can add any such theory description ``key`` into their
rule.

.. code:: {eval-rst}

   .. tip::
     Sometimes, we may want to evaluate a rule provided the perturbative order is within
     a certain range. For example, we may want a rule to be evaluated if the perturbative
     order is strictly less than NLO. This can be done by using directives succeeding the
     :code:`PTO` declaration.

In the above example, one would thus simply use:

.. code:: yaml

     PTO: NLO-

The following are a list of possible directives which can succeed a
``PTO`` declaration: \* ``+`` Evaluate this rule if the theory ``PTO``
is greater than **or equal to** the preceeding PTO \* ``-`` Evaluate
this rule if the theory ``PTO`` is strictly less than the preceeding PTO
\* ``!`` Evaluate this rule if the theory ``PTO`` is not equal to the
preceeding PTO

Examples are:

.. code:: yaml

     PTO: NNLO!
     PTO: N3LO-
     PTO: LO+

If the user doesn’t specify a directive then that implies the rule will
only be evaluated if the declared ``PTO`` matches *exactly* with the
``PTO`` of the theory.

Overwriting filters and default values
--------------------------------------

One can overwrite the default behaviour by adding to the fit runcard.

Custom rules can be added by adding a ``filter_rules:`` namespace in the
fit runcard. This should be a list of rules in the format outlined
above. For example:

.. code:: yaml

   filter_rules:
     - dataset: NMC
       rule: x > 0.2

.. code:: {eval-rst}

   .. warning::
     Adding a :code:`filter_rules` section to the runcard overwrites the default behaviour and does
     **not** append to the default behaviour. This is done intentionally since a rule cannot be 
     overwritten by another rule. By adding the above code snippet, this would be the **only** rule used by
     :code:`vp-setupfit`. As such a bit of copy and pasting may be necessary if one wishes to append a rule.

Similarly the defaults can be overwritten by adding a
``filter_defaults`` namespace to the runcard. For example:

.. code:: yaml

   filter_defaults:
     q2min: 5
     w2min: 10

As in the case of the rules, this overwrites the original defaults and
does not append to them.

.. code:: {eval-rst}

   .. attention::
     To ensure backwards compatibility with old style runcards, if :code:`q2min` and :code:`w2min` are defined
     under the :code:`datacuts` namespace within the runcard, these values are read in and override the default
     values. However, if this overriding occurs, a warning is displayed in standard output.

Examples
--------

Consider the following filter from the ``filters.yaml`` file:

.. code:: yaml

   - dataset: ATLASZPT7TEV
     reason: Avoid the region where resummation effects become important.
     rule: "p_T2 >= 30**2"

this rule applies only to the ``ATLASZPT7TEV`` dataset and keeps all
datapoints with a transverse momentum greater than or equal to 30 MeV.
The reason for the conception of this rule is also provided and we see
that it is due to the fact that datapoints with smaller transverse
momentum will be affected by resummation effects.

Now consider the slightly more complicated example:

.. code:: yaml

   - dataset: CMSDY2D12
     reason: Remove data points for which electroweak corrections are large.
     PTO: NNLO-
     local_variables:
       M: sqrt(M2)
       min_M: 30.0
       max_rapidity: 2.2
     rule: M >= min_M and etay <= max_rapidity

This rule only applies to ``CMSDY2D12``. I wish for the ``rule`` to only
be evaluated provided the ``theory`` perturbative order is **strictly**
less than NNLO (i.e LO or NLO). I check what the process type of
``CMSDY2D12`` is:

.. code:: ipython

   In [1]: from validphys.loader import Loader                                                                                                                                   

   In [2]: l = Loader()                                                                                                                                                          

   In [3]: cd = l.check_commondata("CMSDY2D12")                                                                                                                                  

   In [4]: cd.process_type                                                                                                                                                       
   Out[4]: 'EWK_RAP'

Then cross check this against ``NNPDF.CommonData.kinLabels`` to see that
the relevant kinematic variables are:

::

   'EWK_RAP': ('etay', 'M2', 'sqrts'),

I choose to define custom ``local_variables`` in the form of ``M`` which
is the square root of the invariant mass squared, i.e. just the
invariant mass. Moreover, I define a value for minimum ``M`` and maximum
rapidity which I use in my ``rule`` as cutoff values.

The ``rule`` itself is then self-explanatory, notice however, it is
written in valid ``Python`` syntax. Finally, the reason for the rule is
given which is to cut datapoints which are affected by electroweak
corrections.

As a final example consider the following rule:

.. code:: yaml

   - process_type: DIS_NCP_CH
     reason: |
       Missing higher order corrections to Delta F_IC, the piece that needs
       to be added to the FONLL-C calculation in the case of fitted charm.
     FNS: FONLL-C
     IC: True
     rule: "Q2 > 8"

Instead of this rule applying to one particular dataset, we see it is
applicable to all datasets that have process type ``DIS_NCP_CH``. The
reason for the rule is rather involved and so ``yaml``\ ’s multiline
string syntax is used.

Finally, the user wishes for the ``rule`` to be evaluated **only if**
the theory input has the FONNL-C flavour number scheme and if the theory
uses intrinsic charm. The rule itself is trivial.
