from n3fit.backends import MetaLayer
from n3fit.backends import constraints
from n3fit.backends import operations as op

BASIS_SIZE = 8


class Preprocessing(MetaLayer):
    """
        Applies preprocessing to the PDF.

        This layer generates a factor (1-x)^beta*x^(1-alpha) where both beta and alpha
        are parameters to be fitted.

        Both beta and alpha are initialized uniformly within the ranges allowed in the runcard and
        then they are only allowed to move between those two values (with a hard wall in each side)

        Parameters
        ----------
            `output_dim`: int
                size of the fitbasis
            `flav_info`: list
                list of dicts containing the information about the fitting of the preprocessing
                This corresponds to the `fitting::basis` parameter in the nnpdf runcard.
                The dicts can contain the following fields:
                    `smallx`: range of alpha
                    `largex`: range of beta
                    `trainable`: whether these alpha-beta should be trained during the fit
                                (defaults to true)
            `seed`: int
                seed for the initializer of the random alpha and beta values
    """

    def __init__(
        self,
        output_dim=BASIS_SIZE,
        flav_info=None,
        seed=0,
        initializer="random_uniform",
        **kwargs,
    ):
        self.output_dim = output_dim
        if flav_info is None:
            flav_info = []
        self.flav_info = flav_info
        self.seed = seed
        self.initializer = initializer
        self.kernel = []
        # super(MetaLayer, self).__init__(**kwargs)
        super().__init__(**kwargs)

    def generate_weight(self, weight_name, kind, dictionary):
        """
        Generates weights according to the flavour dictionary and adds them
        to the kernel list of the class

        Parameters
        ----------
            `weight_name`: str
                name to be given to the generated weight
            `kind`: str
                where to find the limits of the weight in the dictionary
            `dictionary`: dict
                dictionary defining the weight, usually one entry from `flav_info`
        """
        limits = dictionary[kind]
        ldo = limits[0]
        lup = limits[1]
        trainable = dictionary.get("trainable", True)
        # Set the initializer and move the seed one up
        initializer = MetaLayer.select_initializer(
            self.initializer, minval=ldo, maxval=lup, seed=self.seed
        )
        self.seed += 1
        # If we are training, constrain the weights to be within the limits
        if trainable:
            weight_constraint = constraints.MinMaxWeight(ldo, lup)
        else:
            weight_constraint = None
        # Generate the new trainable (or not) parameter
        newpar = self.builder_helper(
            name=weight_name,
            kernel_shape=(1,),
            initializer=initializer,
            trainable=trainable,
            constraint=weight_constraint,
        )
        self.kernel.append(newpar)

    def build(self, input_shape):
        # Run through the whole basis
        for flav_dict in self.flav_info:
            flav_name = flav_dict["fl"]
            alpha_name = f"alpha_{flav_name}"
            beta_name = f"beta_{flav_name}"
            self.generate_weight(alpha_name, "smallx", flav_dict)
            self.generate_weight(beta_name, "largex", flav_dict)

        super(Preprocessing, self).build(input_shape)

    def call(self, inputs, **kwargs):
        x = inputs
        pdf_list = []
        for i in range(0, self.output_dim*2, 2):
            pdf_list.append(x ** (1 - self.kernel[i][0]) * (1 - x) ** self.kernel[i+1][0])
        return op.concatenate(pdf_list, axis=-1)
