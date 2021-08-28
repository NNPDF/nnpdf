.. _customplots:

=======================================================
Customizing ``validphys`` plots and other functionality
=======================================================

:code:`validphys` comes with extensive capabilities for producing publication
ready figures and other results. See for example :ref:`pdfplots` or
:ref:`datthcomp`. Here we discuss what to do when these are not quite enough
for a given application, for example when adding extremely specific plots or
producing the final touches for plots for a publication or talk.

Tweaking Matplotlib styles
--------------------------

Some aspects of the appearance of the figure can be customized using
`Matplotlib
stylesheets <https://matplotlib.org/stable/tutorials/introductory/customizing.html>`_.
The default styles we use can be found in the :py:mod:`validphys.mplstyles`
module. A different style file can be passed using the ``--style`` flag of
``validphys`` (and other applications derived from  ``reportengine``). Note
that the styles don't compose and therefore it is advised to copy the default
stylesheet and modify it as needed.

Note that several options haven't been chosen at random. For example the first
few entries in the color cycle are colorblind friendly and such that colors
look OK when stacked with transparency on top of each other, for example for
PDF plots.

Consider adding the functionality to master
-------------------------------------------

New types of plots, as well as stylistic or functional enhancements to existing
code can be :ref:`added to validphys <addvpplots>`, following the :ref:`appropriate
process <rules>`.

This option is strongly recommended and should be the default choice in most
situations. While it requires some initial investment, in coming up with an
appropriate :ref:`design` to make the required feature fit the rest of the
system, as well meeting somewhat :ref:`high coding standards <rules>`, there
are important benefits in exchange: The feature will get a few additional
eyeballs and once merged it will be maintained and kept in sync with the rest
of the code, making runcards using it much more likely to work in the future.
That others can benefit from the work is of course also a good thing.

.. _extramodules:
Hooking ``validphys`` to external code
--------------------------------------

In some situations the requirements for a given plot are rather esoteric and
there is no way to add the functionality to the code economically. In
such cases, external code can be used. Even so, consider upstreaming as much
of the functionality as possible, to get the benefits discussed above.

There are two ways to take advantage of resources produced using the
``validphys`` execution model to process them further.

   * Using the API: It is possible to get some data from ``validphys`` using
     the :ref:`validphys API <vpapi>` and then use it in a script. This affords
     maximum flexibility, as the script can do anything. In exchange runcard
     based input processing or structured output folders aren't readily
     available. Prefer this option for a very small project or when the task
     doesn't fit the execution model of ``validphys`` for some reason.

   * Using extra modules: Additional Python modules or files can be passed to
     ``validphys`` using the ``--extra-modules`` (or ``-x``) flag. The
     functions in these modules then act ``validphys`` providers and can take
     resources from ``validpys`` as input. This approach allows the 
     immediate use of runcards or the default styles. One limitation is that
     there is currently no way of adding production rules or parsers in this
     way. Prefer this for actions that are too difficult to upstream to
     ``validphys``, but should work as if they were internal. A minimal example
     for an external module could be::
         # extra_plots.py

         import matplotlib.pyplot as plt
         from reportengine.figure import figure

         from validphys.commondataparser import load_commondata

         # A simple plot that probably should be in validphys to begin with.

         @figure
         def plot_central_values(commondata):
             fig, ax = plt.subplots()
             ax.plot(load_commondata(commondata).central_values)
             return fig

     The action ``plot_central_values`` can now be used in a runcard:

     .. code:: yaml

               # runcard.py
               dataset_input:
                   dataset: NMC

               actions_:
                   - plot_central_values


    Provided that ``validphys`` is invoked as ``validphys runcard.yaml -x extra_plots.py``.



Note that both of these come at the cost of risking future breakage 
somewhat  as we don't guarantee any sort of stability on the internal
interfaces.

Editing SVG files
-----------------

SVG files store information on figures as sprites and text rather than pixels.
These can then be edited with image editors such as `Inkscape
<https://inkscape.org/>`_. It is possible to edit the text in the figure or
change colors of individual lines.  Note that this is the least maintainable
approach as the modifications need to be applied manually every time the plot
is updated.  However it may be a good way to quickly enhance a plot for a
presentation for example. To produce SVG files, pass  the flag ``--formats
svg`` when invoking ``validphys``.
