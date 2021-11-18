.. _reports:

Generating reports
===================

This section explains how to generate `reports`. These are collections of
plots, tables or other `validphys` outputs which can be uploaded to the
:ref:`vp server <upload>`.

- Reports are implemented as an action of `reportengine`.
- The `report` action takes a `template`
  argument, corresponding to the filename of a template in the
  `Pandoc Markdown format <http://pandoc.org/MANUAL.html#pandocs-markdown>`_, with
  the actions defined with a special syntax discussed below.
- The actions will be resolved as if they were directly specified in the
  configuration file and when all of them are completed, their value
  will be substituted in the template (the `jinja2` library is used for
  the intermediate rendering).

`reportengine` will interpret strings between `{@` and `@}` inside the
templates. There are currently **target** and **with**/**endwith**
tags:

1. Target tags
specify an action to be executed. The possible syntax is:

.. code::

	{@[spec] action_name[(arg1=value, arg2=value)]@}

where `[]` stands for optional syntax. A few conforming examples are:

.. code::

	{@ plot_fancy @}

	{@theory::pdfs plot_fancy@}

	{@plot_fancy(normalize_to=data)@}

The different parts of the specification,
namely mappings, lists of mappings (or special tags implementing that
behaviour) are separated with the `::` operator (resembling the C++
scope resolution operator). Actions will be repeated if the
specification results in multiple namespaces (e.g. one plot per pdf in
the second example above).

2. With/endwith tags
repeat the content between the tags for each namespace in the
specifications. Targets inside the block are repeated and searched for
within each namespace. The syntax of the `with` tag is:

.. code::

	{@with spec@}

and it must be closed by an `endwith` tag

.. code::

  {@endwith@}

Like in the **target** tag, the spec is separated by `::`.

The report action
-----------------

As always, see `validphys --help report` for the most complete
information. The options allow customizing the CSS style or the
template that contains the report itself.

Here we only discuss a couple of interesting flags.

1. The `main` flag

The `main: True` flag can only affect one report per run. It has the
effect of setting the name `index.html`, which comes in handy for
visualizing the uploaded result in the server.

The main flag also tries to open the web browser when the report finishes. The
browser will be chosen according to internal heuristics, by querying system
preferences. These can be overridden by setting the `BROWSER` environment
variable. For example, in text-only environments such as remote clusters, it may
be preferable to just print the URL. This can be achieved by setting the
environment variable to `echo` (for example in the `.bashrc` file):

.. code:: bash

	export BROWSER=echo


2. Displaying math (the `mathjax` flag)

Displaying math on browsers is painful and not without trouble. Pandoc
tries to render the LaTeX math using utf8-characters. This does not
require external dependencies and allows one to work with the text
normally, but is extremely limited (little more than subindexes and
greek letters).

It is possible to set `mathjax: True` to use the
`Mathjax <https://www.mathjax.org/>`_ library. This supports many more
symbols, but is rather slow and requires an external connection in
order to render the math.

Example report template
------------------------

A template that could correspond to the example above is:

.. code::

	NNPDF Report
	============

	{@ description  @}


	PDF plots
	---------

	{@ plot_pdfs @}

	**Normalized**

	{@normalize plot_pdfs  @}


	Train-valid split
	------------------

	{@ plot_training_validation @}

	$\chi^2$
	-------
	{@ with pdfs  @}

	### {@ pdf @}

	{@ experiments_chi2_table @}

	{@ endwith@}

	Experiment plots
	---------------
	{@ with pdfs @}
	###Experiment results for {@pdf@}
	{@with datanorm::experiments@}

	#### {@experiment@}
	{@experiment plot_fancy @}
	{@ endwith @}
	{@ endwith @}


First we are writing a verbatim Markdown title. Next we are asking for
a variable named "`description`" to be computed and later substituted
right below (it is obtained from the fit config file, as seen in the
template). Then we are computing absolute and normalized PDF plots
(`normalize` is an arbitrary string that is defined in the config file
to normalize to the first PDF). We then plot the training and
validation :math:`\chi^2` of each replica in the fit. Next we compute the
:math:`\chi^2` for each experiment, and produce a separate table and heading
for each PDF in `pdfs` (note that LaTeX math syntax is allowed).
Finally we produce, for each pdf and for each experiment, a set of
data-theory comparison plots (which in turn are repeated for each
dataset in the experiment).

Customizing how things look in the report
-----------------------------------------

By default, the `str()` method will be applied to objects that appear
in the report. If you want a custom behaviour, declare
a custom `as_markdown` property for your objects. It should return
a string in Pandoc Markdown describing your object. Raw HTML is
also allowed (although that decreases the compatibility, e.g. if we
decide to output LaTeX instead of HTML in the future).
