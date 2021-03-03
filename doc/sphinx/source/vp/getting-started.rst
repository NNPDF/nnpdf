Getting started with validphys
==============================

To use ``validphys`` you must provide an input runcard which includes

* The resources you need (PDFs, fits, etc.)
* The actions (functions) you would like to be carried out
* Additional flags and parameters for fine-tuning
* Metadata describing the author, title and keywords

To get an idea of the layout, :ref:`vpexamples` details the example runcards that can be found in 
`this folder <https://github.com/NNPDF/nnpdf/tree/master/validphys2/examples>`_. The :ref:`tutorials`
section also takes you through how to make runcards for various tasks.

Once you have created a runcard (e.g. ``runcard.yaml``), simply run

.. code:: 

   validphys runcard.yaml

to set the ball rolling.

Another useful command to be aware of is ``vp-comparefits - i``, which launches an interactive
session to compare two fits. See the tutorial :ref:`compare-fits` for more information.

For more tailored analysis, the API provides a high level interface to the code, allowing you to 
extract objects and play around with them. See :ref:`vpapi`.

Finally, the ``validphys --help`` command can give you information on modules and specific actions, e.g.

.. code:: 
   
   	$ validphys --help fits_chi2_table
   
   	fits_chi2_table

	Defined in: validphys.results

	Generates: table

	fits_chi2_table(fits_total_chi2_data, fits_datasets_chi2_table,
	  fits_groups_chi2_table, show_total: bool = False)

	Show the chi² of each and number of points of each dataset and
	experiment of each fit, where experiment is a group of datasets
	according to the `experiment` key in the PLOTTING info file, computed
	with the theory corresponding to the fit. Dataset that are not
	included in some fit appear as `NaN`



	The following additionl arguments can be used to control the
	behaviour. They are set by default to sensible values:

	  show_total(bool) = False
	  per_point_data(bool) = True [Used by fits_groups_chi2_table]
	  
	  


