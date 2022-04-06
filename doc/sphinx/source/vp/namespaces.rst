.. _namespaces:

Namespaces
==========

In this section we explain a little bit about namespaces. Namespaces are used frequently in
`validphys` runcards because these runcards are based on `reportengine`, for which namespaces
are a central concept.

A namespace is a stack of python dictionaries which are indexed by a tuple called the
`namespace specification (nsspec)`. A user gives some inputs in terms of `fuzzyspecs` and
these are translated into `nnspecs`. This is mostly an advanced internal implementation detail,
but it is important in order to understand how several features work. Also, the abstraction leaks
into user-facing features such as the :ref:`collect function <collect>`.

Namespace specifications
------------------------

An nsspec is a tuple of an arbitrary number of elements. Each element
in the tuple corresponds to one extra stack layer in depth (*"stack
frame"*). The elements of the tuple can be either:

 - Names of mappings.
 - Names of objects that have an `as_namespace` method
 - Tuples of the form (name of list of mappings, index).

The scope rules are similar to those of C: the lookup of a value is
done first looking at the inner frame and then at the outer ones,
until a match is found.

Consider the following example:

.. code:: yaml

	First:
	   pdf: NNPDF40_nlo_as_01180
	   normalize_to: None
	   use_cuts: "nocuts"

	Second:
	   pdf: NNPDF40_nnlo_as_01180
	   normalize_to: NNPDF40_nnlo_as_01180

	cutspecs:
	 - {use_cuts: "nocuts"}
 	 - {use_cuts: "fromfit"}


Given the input above, we could form the following `nsspec`:

.. code:: python

	('Second', ('cutspecs', 0))

This would correspond to a namespace where we have the following
symbols available:

- `use_cuts` (set to `"nocuts"`) from `cutspecs`.
- `pdf` and `normalize_to` (set to NNLO) from `second`.
- `First`, `Second` and `cutspecs` from the root namespace.

We could also form the specification:

.. code:: python

	(('cutspecs', 1), 'First')

Because the innermost specification is last, the value of `use_cuts`
is `"nocuts"`.

The function `reportengine.namespaces.resolve(ns, nsspec)` returns
a mapping (in particular it is a modified version of
`collections.ChainMap`) that implements exactly this behaviour. It is
used extensively thorough `reportengine`.

.. _fuzzyspecs:

Fuzzyspecs
----------

The namespace specifications as described above is not what
the user typically enters. Instead, the typical user input is what in
the code is labelled *fuzzyspec*. A fuzzyspec is like an nsspec except
that the lists of mappings are entered by name and not by a tuple
(name, index). A fuzzyspec resolves to one or more nsspecs. For
example, given the fuzzyspec:

.. code:: python

	('Second', 'cutspecs')

and the input above, it gets expanded into two nsspecs:

.. code:: python

	('Second', ('cutspecs', 0))
	('Second', ('cutspecs', 1))

corresponding to each of the two mappings in cutspecs.

The `as_namespace` method
-------------------------

An object can customize how it is viewed as a reportengine namespace.
This is done by defining a method called `as_namespace`, which takes no
arguments and should return either a mapping or a list of mappings.
This is used to implement parsing lists automatically.
