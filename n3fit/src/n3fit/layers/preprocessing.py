from typing import List, Optional

from n3fit.backends import MetaLayer, MultiInitializer, constraints
from n3fit.backends import operations as op


class Preprocessing(MetaLayer):
    """
    Computes preprocessing factor for the PDF.

    This layer generates a factor (1-x)^beta*x^(1-alpha) where both beta and alpha
    are model paramters that can be trained. If feature scaling is used, the preprocessing
    factor is x^(1-alpha).

    Alpha is initialized uniformly within the ranges allowed in the runcard and
    then it is only allowed to move between those two values (with a hard wall in each side)

    Alpha and, unless feature scaling is used, beta are initialized uniformly within
    the ranges allowed in the runcard and then they are only allowed to move between those two
    values (with a hard wall in each side)

    Parameters
    ----------
        replica_seeds: List[int]
            list of pre replica seeds for the initializer of the random alpha and beta values
        flav_info: list
            list of dicts containing the information about the fitting of the preprocessing factor
            This corresponds to the `fitting::basis` parameter in the nnpdf runcard.
            The dicts can contain the following fields:
                `smallx`: range of alpha
                `largex`: range of beta
                `trainable`: whether these alpha-beta should be trained during the fit
                            (defaults to true)
        large_x: bool
            Whether large x preprocessing factor should be active
    """

    def __init__(
        self,
        replica_seeds: Optional[List[int]],
        flav_info: Optional[list] = None,
        large_x: bool = True,
        **kwargs,
    ):
        if flav_info is None:
            raise ValueError(
                "Trying to instantiate a preprocessing factor with no basis information"
            )
        self.flav_info = flav_info
        self.replica_seeds = replica_seeds
        self.large_x = large_x
        self.num_replicas = len(replica_seeds)

        self.alphas = []
        self.betas = []
        super().__init__(**kwargs)

    def generate_weight(self, name: str, kind: str, dictionary: dict, set_to_zero: bool = False):
        """
        Generates weights according to the flavour dictionary

        Parameters
        ----------
            name: str
                name to be given to the generated weight
            kind: str
                where to find the limits of the weight in the dictionary
            dictionary: dict
                dictionary defining the weight, usually one entry from `flav_info`
            set_to_zero: bool
                set the weight to constant 0
        """
        constraint = None
        if set_to_zero:
            single_replica_initializer = MetaLayer.init_constant(0.0)
            trainable = False
        else:
            minval, maxval = dictionary[kind]
            trainable = dictionary.get("trainable", True)
            # Seeds will be set in the wrapper MultiInitializer
            # Note: keras 3 interprets a seed of 0 as None and replaces it with a random seed,
            # so we set it to 1 here and subtract it later when the replica seed is set
            single_replica_initializer = MetaLayer.select_initializer(
                "random_uniform", minval=minval, maxval=maxval
            )
            # If we are training, constrain the weights to be within the limits
            if trainable:
                constraint = constraints.MinMaxWeight(minval, maxval)

        initializer = MultiInitializer(single_replica_initializer, self.replica_seeds, base_seed=0)
        # increment seeds for the next coefficient
        self.replica_seeds = [seed + 1 for seed in self.replica_seeds]

        # Generate the new trainable (or not) parameter
        newpar = self.builder_helper(
            name=name,
            kernel_shape=(self.num_replicas, 1),
            initializer=initializer,
            trainable=trainable,
            constraint=constraint,
        )
        return newpar

    def build(self, input_shape):
        # Run through the whole basis
        for flav_dict in self.flav_info:
            flav_name = flav_dict["fl"]
            alpha_name = f"alpha_{flav_name}"
            self.alphas.append(self.generate_weight(alpha_name, "smallx", flav_dict))
            beta_name = f"beta_{flav_name}"
            self.betas.append(
                self.generate_weight(beta_name, "largex", flav_dict, set_to_zero=not self.large_x)
            )

        super(Preprocessing, self).build(input_shape)

    def call(self, x):
        """
        Compute preprocessing prefactor.

        Parameters
        ----------
            x: tensor(shape=[1,N,1])

        Returns
        -------
            prefactor: tensor(shape=[1,R,N,F])
        """
        # weight tensors of shape (R, 1, F)
        alphas = op.stack(self.alphas, axis=-1)
        betas = op.stack(self.betas, axis=-1)

        x = op.batchit(x, batch_dimension=0)

        return x ** (1 - alphas) * (1 - x) ** betas
