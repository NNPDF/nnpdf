.. _releases:
Releases and compatibility policy
=================================

We follow a rolling development model where the tip of the `master branch
<https://github.com/NNPDF/nnpdf/tree/master>`_ is expected to be stable, tested
and correct. Binary packages for the latest commit on the branch, with
appropriate version information are :ref:`generated automatically<CI>` and can
be :ref:`readily installed<conda>`. In general the version of the code should be
preferred for producing new results, but see the compatibility policy below. The
main results, such as NNPDF 4.0 :cite:p:`nnpdf40` will be produced with a frozen
:ref:`tag <tags>`, a :ref:`conda environment <conda>` and a :ref:`docker image
<docker>` so that they can be :ref:`reproduced <reproduce40>` entirely.

.. _tags:
Tags
----

The code is tagged to contextualize the versioning, mark significant
developments and to mark versions used to produce main results. The
significant releases since the code was made public are:

`Version 4.0.7 <https://github.com/NNPDF/nnpdf/releases/tag/4.0.7>`_
    A pre-release for the QED, MHOU and N3LO papers
`Version 4.0.6 <https://github.com/NNPDF/nnpdf/releases/tag/4.0.6>`_
    The last version that uses C++ objects in n3fit and validphys, and apfel for
    PDF evolution.
`Version 4.0.5 <https://github.com/NNPDF/nnpdf/releases/tag/4.0.5>`_
    The last version to support legacy genetic algorithms fits based on C++.
`Version 4.0.4 <https://github.com/NNPDF/nnpdf/releases/tag/4.0.4>`_
    Incremental bugfix and enhancementent release: most importantly fixing an
    error in the definition of integrated luminosity (`Github issue #1442
    <https://github.com/NNPDF/nnpdf/issues/1442>`_), as well as adding
    ancillary features for the NNPDF4.0 analysis and results, and extending
    the documentation.
`Version 4.0.3 <https://github.com/NNPDF/nnpdf/releases/tag/4.0.3>`_
    Used to :ref:`produce <reproduce40>` the NNPDF 4.0 :cite:p:`nnpdf40`
    fits.

.. _compatibility_policy:
Compatibility Policy
--------------------

Fit runcards, and Physics results
````````````````````````````````````

We attempt to maintain strict compatibility for all published results for at
least one major PDF release cycle, in order to be able to reproduce the
current published release (while the new one is under development) and compare
new and old results. For example the code that produces NNPDF 4.0 should be
able to reproduce the results for 3.1 as well. Once NNPDF 4.0 is released, new
developments in the code are allowed to break compatibility with 3.1, but
should maintain it with 4.0 (until after 4.1 would be released).

The baseline expectation for fits is that a
given runcard for a published PDF set is able to produce an equivalent PDF. If
bugs are discovered in experimental data or theory predictions, new fixed
copies of the data would be made with different names, while keeping the old
ones in such a way that they are selected with the original runcards. Minor
differences in the result can happen as a consequence of small enhancements or
differences in random number generation, although that would be avoided when
practicable. If significative changes are necessary, for example as a result of
discovering a bug, this will be clearly indicated by a tagged release.


Analysis runcards used for published results are expected to be able to produce
the *same physics*, while bugfixes (that don't affect fits) or changes in
presentation can happen in between. Similarly, important enough bugfixes will
be marked by a tag.

Internal interfaces
`````````````````````

We follow a `"Linux Kernel"
<https://en.wikipedia.org/wiki/Linux_kernel_interfaces#In-kernel_APIs>`_
approach to internal interfaces, which do not affect the content of runcards.
This means that there is no expectation of stability at all and these
interfaces can change arbitrarily at every commit without any particular
notice. If you wish that code such as :ref:`extra modules<extramodules>` is
maintained and kept in working order with newer updates, it is highly
suggested to :ref:`contribute it to the main repository <rules>`,
along with appropriate tests and documentation. Otherwise you are on your
own.
