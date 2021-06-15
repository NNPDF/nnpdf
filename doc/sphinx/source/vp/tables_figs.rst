Producing tables and figures
============================

Producing figures
-----------------

In order to produce figures, just decorate your functions returning
`matplotlib` `Figure` objects  with the `reportengine.figure.figure`
function, e.g.:

.. code:: python

	@figure
	def plot_p_alpha(p_alpha_study):
	   fig, ax = plt.subplots()
	   #Plot something
	   ...
	   return fig

This will take care of the following:

 - Saving the figures with a nice, unique name to the output folder,
   in the formats specified by the user.

 - Closing the figures to save memory.

 - Making sure figures are properly displayed in reports.

There is also the `figuregen` decorator for providers that are
implemented as generators that yield several figures (see e.g. the
implementation of `plot_fancy`). Apart from just the figure, they yield
a tuple (prefix, figure) where the prefix will be used in the
filename.

Producing tables
----------------

These work similarly to producing figures, as described
above. Instead use the `@table` and `@tablegen` decorators.

Tables will be saved in the CSV formats.
