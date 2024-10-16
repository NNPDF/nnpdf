.. _tut_report:

How to generate a report
========================

Suppose that we want to generate a custom report that includes plots and
statistics that are not included as part of the report generated by
:ref:`vp-comparefits <compare-fits>`. We may be lucky enough to find an example
runcard that produces what we need in
`validphys2/examples <https://github.com/NNPDF/nnpdf/tree/master/validphys2/examples>`_.
However, we may need to write our own own `yaml <https://yaml.org/>`_ runcard
from scratch.

Suppose we want to have histograms of the χ2 per replica for some set of
experimental datasets.

The calling of :ref:`validphys <vp-index>` actions are used as normal. The action we
are looking for is
`plot_chi2dist <https://github.com/NNPDF/nnpdf/blob/d79059975e4ef97063c6bdd9f19dfb908586e453/validphys2/src/validphys/dataplots.py#L50>`_.
Here's an example report that does what we're looking for:

.. code:: yaml

  meta:
    title: Distribution plots per replica across experiments
    author: Shayan Iranipour
    keywords: [chi2, replica, distribution, DISonly]

  fit: NNPDF31_nnlo_as_0118_DISonly

  pdf:
    from_: "fit"

  experiments:
    from_: "fit"

  theoryid: 53

  use_cuts: "fromfit"

  template_text: |
    # Histograms of χ2
    ## DIS only distributions
    {@experiments::experiment plot_chi2dist@}

  actions_:
    - report(main=True)

The ``report(main=True)`` command is what generates the report. We can customize
the formatting of the report using
`markdown <https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet>`_
syntax. Note for example that ``# Histograms of χ2`` gives an appropriate title
to the section where we will find our plot.

If the ``template_text`` section of the runcard becomes large and unwieldy, it may
be preferable to put the information from this section in a separate file. In
such cases one can create a ``markdown`` template file, usually called ``template.md``
such as

.. code:: md

  # Histograms of χ2
  ## DIS only distributions
  {@experiments::experiment plot_chi2dist@}

and change the `template_text` section of the runcard to the following

.. code:: yaml

  template: template.md

where this assumes that ``template.md`` is in the same folder as that in which you
execute the ``validphys`` command.

The ``meta`` field is important for retrieving the report once it has been
uploaded to the `validphys server <https://vp.nnpdf.science/>`_. ``title`` and
``author`` are fields that appear when browsing through reports while the
``keywords`` allow quick retrieval of the report by the search functionality of
the server. Setting appropriate keywords is especially important when working on
a big project, within which it is likely that many reports will be produced. In
such cases a ``keyword`` should be chosen for the project and set in each uploaded
report.
