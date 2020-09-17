Backend selection
=================

One of the advantages of ``n3fit`` is the possibility to select different backends
for training the neural network.

Currently implemented in ``n3fit`` we can find:

- [Keras-Tensorflow](#keras-tensorflow)
- [Evolutionary Algorithms](#evolutionary-algorithms)

Keras-Tensorflow
----------------
The main and default backend for ``n3fit`` is the [TensorFlow](https://www.tensorflow.org/)
library developed by Google.


Evolutionary Algorithms
-----------------------
As part of the validation process followed during the development of ``n3fit`` we developed the
[evolutionary-keras](https://evolutionary-keras.readthedocs.io) [[1](https://arxiv.org/abs/2002.06587)] library.

This library extends [Tensorflow](#keras-tensorflow) for its use with evolutionary algorithms such as the [NGA](https://evolutionary-keras.readthedocs.io/en/latest/optimizers.html#nodal-genetic-algorithm-nga)
used in previous versions of [NNPDF](http://arxiv.org/abs/1410.8849) or the [CMA-ES](https://evolutionary-keras.readthedocs.io/en/latest/optimizers.html#covariance-matrix-adaptation-evolution-strategy-cma-es) used in some [spin-off studies](https://arxiv.org/abs/1706.07049).

The evolutionary strategies can be used by adding the keyword ``backend: evolutionary_keras`` to the ``fitting`` section in the runcard.

```yaml
fitting:
  backend: evolutionary_keras
  optimizer:
    optimizer_name: 'NGA'
    sigma_init: 15
    population_size: 80
    mutation_rate: 0.05
```
