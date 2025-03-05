.. _releases:
Releases and compatibility policy
=================================

Development occur in the tip of the `master branch <https://github.com/NNPDF/nnpdf/tree/master>`_
while we aim for this branch to be stable, tested and correct, this is not guaranteed.
Premade packages are available for the latest tag
:ref:`generated automatically<CI>` and can be :ref:`readily installed<conda>`.
See the compatibility policy below. The
main results, such as NNPDF 4.0 :cite:p:`nnpdf40` will be produced with a frozen
:ref:`tag <tags>`, a :ref:`conda environment <conda>` and a :ref:`docker image
<docker>` so that they can be reproduced entirely.

.. _tags:
Tags
----

The code is tagged to contextualize the versioning, mark significant
developments and to mark versions used to produce main results. The
significant releases since the code was made public are:

`Version 4.0.9 <https://github.com/NNPDF/nnpdf/releases/tag/4.0.9>`_
    Release for 4.0 `N3LO <https://arxiv.org/abs/2402.18635>`_;
    last release fully backwards-compatible with 4.0 pipeline. 4.0 runcards will still work but
    external tools, and data and theory not used in the 4.0 family of fits will no longer be
    guaranteed to work from 4.0.10 onwards Last release compatible with the old commondata format
    and that accepts apfel as evolution code.
`Version 4.0.8 <https://github.com/NNPDF/nnpdf/releases/tag/4.0.8>`_
    Release for the `QED <https://arxiv.org/abs/2401.08749>`_ and `MHOU <https://arxiv.org/abs/2401.10319>`_ papers.
    Miscellaneous bugfixes and small QOL improvements. See the whole list of changes in the release description.
    This is the last version that uses cmake for installation.
`Version 4.0.7 <https://github.com/NNPDF/nnpdf/releases/tag/4.0.7>`_
    Intermediate release with miscelanous improvements in preparation for 4.0.8.
    Development release with experimental features
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
    Used to produce the NNPDF 4.0 :cite:p:`nnpdf40` fits.

.. _compatibility_policy:
Compatibility Policy
--------------------

Fit runcards, and Physics results
`````````````````````````````````

We attempt to maintain strict compatibility for all published results for at
least one major PDF release cycle, in order to be able to reproduce the
current published release (while the new one is under development) and compare
new and old results. For example the code that produces NNPDF 4.0 should be
able to reproduce the results for 3.1 as well. Once NNPDF 4.0 is released, new
developments in the code are allowed to break compatibility with 3.1, but
should maintain it with 4.0 (until the 4.1 cycle and so on).

The baseline expectation for fits is that a given runcard for a published PDF set
is able to produce an equivalent PDF when using the same tag of the code.
Due to storage limitations, the theory object (FastKernel tables and EKOs)
are not versioned, therefore if bugs are discovered in theory predictions
new fixed theories will be made available which might change some of these results.

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
along with appropriate tests and documentation.:
