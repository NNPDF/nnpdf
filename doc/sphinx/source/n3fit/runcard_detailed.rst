================================
``n3fit`` runcard detailed guide
================================

In this section we fine-grain the explanation of the different parameters that enter the runcard.

- :ref:`preprocessing-label`
- :ref:`trval-label`
- :ref:`networkarch-label`
- :ref:`optimizer-label`
- :ref:`positivity-label`
- :ref:`otheroptions-label`


.. _preprocessing-label:

Preprocessing
-------------
The behaviour of the preprocessing in the ``n3fit`` code is controlled, as in the old ``nnfit`` code, through the ``fitting:basis`` parameter of the nnpdf runcard.

The preprocessing factor applied to every flavour of the basis is:

.. math::

    x ^ {1 - \alpha} (1 - x) ^{\beta}


This parameter accepts a list of the size of the chosen basis with a number of parameter for each flavour. The parameters  used in ``n3fit`` are:

- ``fl``: name of the flavour, this name will be use to define the name of the weights as ``alpha_{fl}`` and ``beta_{fl}``.
- ``smallx``: range of the ``alpha``
- ``largex``: range of the ``beta``
- ``trainable``: sets the flavour basis to be trainable or not, defaults to ``True``

Setting the ``trainable`` flag to ``False`` is equivalent to recovering the old behaviour of ``nnfit``.

.. code-block:: yaml

    fitting:
        basis:
            # smallx, largex: preprocessing ranges
            - { fl: sng, smallx: [1.05,1.19], largex: [1.47,2.70], trainable: False }
            - { fl: g,   smallx: [0.94,1.25], largex: [0.11,5.87], trainable: False }
            - { fl: v,   smallx: [0.54,0.75], largex: [1.15,2.76], trainable: False }
            - { fl: v3,  smallx: [0.21,0.57], largex: [1.35,3.08] }
            - { fl: v8,  smallx: [0.52,0.76], largex: [0.77,3.56], trainable: True }
            - { fl: t3,  smallx: [-0.37,1.52], largex: [1.74,3.39] }
            - { fl: t8,  smallx: [0.56,1.29], largex: [1.45,3.03] }
            - { fl: cp,  smallx: [0.12,1.19], largex: [1.83,6.70] }


.. _trval-label:

Training / Validation split
---------------------------
The fraction of events that are considered for the training and validation sets is defined by the ``frac`` key in the ``experiment:dataset`` parameter of the nnpdf runcard. A fraction of ``X`` means that ``X`` of the event will go into the training set while ``1-X`` will enter the validation set for that dataset.

.. code-block:: yaml

    experiments:
    - experiment: ALL
        datasets:
        - { dataset: SLACP, frac: 0.8}
        - { dataset: NMCPD, frac: 0.8 }      
        - { dataset: CMSJETS11,     frac: 0.8, sys: 10 }

It is possible to run a fit with no validation set by setting the fraction to ``1.0``, in this case the training set will be used as validation set.


.. _networkarch-label:

Network Architecture
--------------------
There are different network architectures implemented in ``n3fit``.
Which can be selected by changing the ``fitting:parameters::layer_type`` parameter in the runcard.
All layer types implement the ``nodes_per_layer``, ``activation_per_layer`` and ``initializer`` parameters.

.. code-block:: yaml

    fitting:
        parameters:
            nodes_per_layer: [5, 3, 8]
            activation_per_layer: ['tanh', 'tanh', 'linear']
            layer_type: 'dense_per_flavour'
            initializer: 'glorot_normal'

- **One single network** (``layer_type: dense``):

  Extra accepted parameters:
    - `dropout`: float
        see `keras dropout <https://keras.io/layers/core/#dropout>`_
    - `regularizer`: str
        see `keras regularizers <https://keras.io/regularizers/>`_
    - `regularizer_args`: dict
        choice arguments for the `regularizer`

In this mode all nodes are connected with all nodes of the next layer. In this case there is one single network which take as input the value of ``x`` (and ``log(x)``) and outputs all different flavours.

In this case the ``nodes_per_layer`` parameter represents the nodes each one of these layers has. For instance ``[40, 20, 8]`` corresponds to a network where the first layer is a matrix ``(2x40)`` (the input is ``x, log(x)``), the second layer is a matrix ``(40x20)`` and the third and final one ``(20x8)``.

- **One network per flavour** (``layer_type: dense_per_flavour``):

This mode is designed to behave as the methodology for NNPDF before 3.1 where each flavour has a separated identical network. 

In this case the ``nodes_per_layer`` parameter represents the nodes each layer of each flavour has. For instance ``[5, 3, 8]`` means that the first step is a list of 8 layers of shape ``(2x5)``, while the second layer is again a list that matches the previous one (i.e., 8 layers) with layers of shape ``(5x3)`` while the last layer has two task. The output of each layer should be one single element (i.e., 8 ``(3x1)`` layers) and then concatenate them all so that the final output of the neural network will be a 8-elements tensor. A report comparing the ``dense`` and ``dense_per_flavour`` architectures can be found  `here <https://vp.nnpdf.science/q6Rm1Q_rTguJwKsLOZFoig==/>`_


.. _optimizer-label:

Optimizer
---------

One of the most important parameters defining the training of the Neural Network is the choice
of optimizer (and its corresponding options).

.. code-block:: yaml

    fitting:
        parameters:
            optimizer:
              optimizer_name: 'Adadelta'
              learning_rate: 1.0
              clipnorm: 1.0


The full list of optimizers accepted by the ``n3fit`` and their arguments
can be checked in the `MetaModel <https://github.com/NNPDF/nnpdf/blob/master/n3fit/src/n3fit/backends/keras_backend/MetaModel.py>`_ file.



.. _positivity-label:

Positivity
----------

In ``n3fit`` the behavior of the positivity observables has changed with respect to ``nnfit``.
In ``nnfit`` the loss due to the positivity observable was multiplied by a ``poslambda`` for each observable, defined in the runcard as:

.. code-block:: yaml

    positivity:
      posdatasets:
        - {dataset: POSF2U, poslambda: 1e6}


This behavior was found to be very inefficient for gradient descent based strategies and was exchanged for a dynamical Lagrange multiplier.
The dynamical multiplier is defined in terms of a initial value and a multiplier to be applied each 100 epochs.
Both the initial value and the 100 epochs multiplier are defined as an optional ``positivity`` dictionary alongside the hyperparameters of
the Neural Network as:

.. code-block:: yaml

    fitting:
        parameters:
            positivity:
              multiplier: 1.05
              initial: 14.5
              
Note that by defining the positivity in this way all datasets will share the same Lagrange multiplier.

It is also possible to not define the positivity hyperparameters (or define them only partially).
In this case ``n3fit`` will set the initial Lagrange multiplier as ``initial`` (default: 1.0)
while the ``multiplier`` will be such that after the last epoch the final Lagrange multiplier 
equals the ``poslambda`` defined for the dataset.

Note that the effect of the positivity constraints on the fit are still being studied and so this cannot be considered final.
     
              
.. _otheroptions-label:

Other options
-------------


.. code-block:: yaml

    fitting:
        parameters:
            threshold_chi2: 4.0

- ``threshold_chi2``: sets a maximum validation :math:`\chi2` for the stopping to activate. Avoids (too) early stopping.

