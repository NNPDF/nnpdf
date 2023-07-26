from n3fit.backends import MetaLayer, constraints
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
        seed: int
            seed for the initializer of the random alpha and beta values
    """

    def __init__(
        self,
        flav_info: list = None,
        seed: int = 0,
        large_x: bool = True,
        **kwargs,
    ):
        if flav_info is None:
            raise ValueError(
                "Trying to instantiate a preprocessing factor with no basis information"
            )
        self.flav_info = flav_info
        self.seed = seed
        self.initializer = "random_uniform"
        self.large_x = large_x
        self.alphas = []
        self.betas = []
        super().__init__(**kwargs)

    def generate_weight(self,
        name: str,
        kind: str,
        dictionary: dict,
        set_to_zero: bool = False
        ):
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
            initializer = MetaLayer.init_constant(0.0)
            trainable = False
        else:
            minval, maxval = dictionary[kind]
            trainable = dictionary.get("trainable", True)
            # Set the initializer and move the seed one up
            initializer = MetaLayer.select_initializer(
                self.initializer, minval=minval, maxval=maxval, seed=self.seed
            )
            self.seed += 1
            # If we are training, constrain the weights to be within the limits
            if trainable:
                constraint = constraints.MinMaxWeight(minval, maxval)

        # Generate the new trainable (or not) parameter
        newpar = self.builder_helper(
            name=name,
            kernel_shape=(1,),
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
            self.betas.append(self.generate_weight(beta_name, "largex", flav_dict, set_to_zero=not self.large_x))

        super(Preprocessing, self).build(input_shape)

    def call(self, x):
        """
        Compute preprocessing prefactor.

        Parameters
        ----------
            x: tensor(shape=[1,N,1])

        Returns
        -------
            prefactor: tensor(shape=[1,N,F])
        """
        alphas = op.stack(self.alphas, axis=1)
        betas = op.stack(self.betas, axis=1)

        return x ** (1 - alphas) * (1 - x) ** betas
