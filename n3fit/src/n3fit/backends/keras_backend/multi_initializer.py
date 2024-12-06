"""
    Extend the ``Initializer`` from Keras to initialize an arbitrary number of replicas,
    with a behaviour equal to initiailizing a bunch of replicas and then stacking them.
"""

from keras.initializers import Initializer

from .operations import stack


class MultiInitializer(Initializer):
    """
    Multi replica initializer that exactly replicates a stack of single replica initializers.

    Weights are stacked on the first axis, and per replica seeds are added to a base seed of the
    given single replica initializer.

    Parameters
    ----------
        single_initializer: Initializer
            Initializer class for the kernel.
        replica_seeds: List[int]
            List of seeds per replica for the kernel initializer.
        base_seed: int
            Base seed for the single replica initializer to which the replica seeds are added.
    """

    def __init__(self, single_initializer: Initializer, replica_seeds: list[int], base_seed: int):
        self.initializer_class = type(single_initializer)
        self.initializer_config = single_initializer.get_config()
        self.base_seed = base_seed
        self.replica_seeds = replica_seeds

    def __call__(self, shape, dtype=None, **kwargs):
        shape = shape[1:]  # Remove the replica axis from the shape.
        per_replica_weights = []
        for replica_seed in self.replica_seeds:
            if "seed" in self.initializer_config:
                self.initializer_config["seed"] = int(self.base_seed + replica_seed)
            single_initializer = self.initializer_class.from_config(self.initializer_config)

            per_replica_weights.append(single_initializer(shape, dtype, **kwargs))

        return stack(per_replica_weights, axis=0)
