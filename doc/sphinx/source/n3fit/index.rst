.. _n3fitindex:

Fitting code: ``n3fit``
=======================

-  ``n3fit`` is the next generation fitting code for NNPDF developed by the
   N3PDF team :cite:p:`Carrazza:2019mzf`
-  ``n3fit`` is responsible for fitting PDFs from NNPDF4.0 onwards.
-  The code is implemented in python using `Keras <https://keras.io/>`_ and can run with `Tensorflow <https://www.tensorflow.org>`_ (default) or `pytorch <https://pytorch.org>`_ (with the environment variable ``KERAS_BACKEND=torch``).
-  The sections below are an overview of the ``n3fit`` design.


.. important::
    If you just want to know how to run a fit using ``n3fit``, head to :ref:`n3fit-usage`.


.. toctree::
   :maxdepth: 2

   methodology
   runcard_detailed
   hyperopt
