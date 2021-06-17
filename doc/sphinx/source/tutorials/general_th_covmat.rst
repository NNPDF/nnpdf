Including a general theory covariance matrix in a fit
=====================================================
This tutorial explains how to include an externally constructed theory covariance 
matrix (theory covmat) in a fit. 

.. warning::
   Theory covariance matrices are currently only supported in the legacy fitting code,
   ``nnfit``, which was used to produce NNPDF3.1 fits. See :ref:`nnfit-usage`.

.. note::
   Scale variation (MHOU) covariance matrices are already implemented in ``validphys``
   in the `theorycovariance <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/theorycovariance/>`_ module. 
   See the tutorial on how to include them.
   
Sometimes we would like to see the impact of a certain theory covariance matrix 
on the fit. For example

-  Higher twist uncertainties
-  Top mass uncertainties
-  Nuclear or deuteron uncertainties

.. note::
    Currently nuclear and deuteron uncertainties are implemented via 
    `buildmaster <https://github.com/NNPDF/nnpdf/tree/master/buildmaster/>`_,
    using additional systematics.
    
Adding a user covmat can be done using the ``user_covmat`` action in
`theorycovariance/construction.py <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/theorycovariance/construction.py>`_.

Instructions
------------
1. Save your covmat as a .csv file of a pandas DataFrame Multiindexed 
   in the same way as data in ``validphys``. See the actions ``groups_covmat`` and 
   ``groups_index``
   in `results.py <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/results.py>`_ as an example. You can use 
   `pandas.DataFrame.to_csv <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_csv.html>`_ to do this.

2. Put the covmat in a folder and :ref:`vp-upload <vpupload>` it. 

.. warning:: 
    Make a note of the upload address returned to you, but without the initial
    part of the address, i.e. you should save
    "https://vp.nnpdf.science/IeGM9CY8RxGcb5r6bIEYlQ==/shrek_covmat.csv" 
    as "IeGM9CY8RxGcb5r6bIEYlQ==/shrek_covmat.csv"

3. In the runcard under ``theorycovmatconfig`` you need to add the 
   following (using the address above as an example)

.. code:: yaml

    ############################################################################
    theorycovmatconfig:
      use_user_uncertainties: True
      user_covmat_path: "IeGM9CY8RxGcb5r6bIEYlQ==/shrek_covmat.csv"
      use_thcovmat_in_sampling: True
      use_thcovmat_in_fitting: True		
    ############################################################################

The flags ``use_thcovmat_in_fitting`` and ``use_thcovmat_in_sampling`` specify
where to use the theory covmat in the code. There are two possible places:
the fitting (i.e. \\(\\chi^2\\) minimiser) and the sampling (i.e. pseudodata
generation). The default is ``True`` for both.

.. warning::
      Changing either of these to ``False`` will affect the fit outcome and should
      be avoided unless you know what you are doing.
      
4. Make sure that datasets are grouped under one big experiment called "BIGEXP", 
   just like in :ref:`vptheorycov-index`.
   
5. For an example runcard, see `here <https://github.com/NNPDF/nnpdf/tree/master/validphys2/examples/fit_with_user_covmat.yaml.>`_

.. warning::
	Theory covariance matrices can currently only be included using the legacy `nnfit` code. When 	running `vp-setupfit` you need to include the `--legacy` flag

Including both scale variation uncertainties and user uncertainties
-------------------------------------------------------------------
User uncertainties and scale variation uncertainties are included independently.
By default neither are included. To include both
see the separate tutorial on scale variation uncertainties and use the 
union of the contributions in ``theorycovmatconfig``.	For an example runcard see `here <https://github.com/NNPDF/nnpdf/tree/master/validphys2/examples/fit_with_sv_and_user_covmat.yaml.>`_
