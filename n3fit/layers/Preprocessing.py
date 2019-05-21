from backends import MetaLayer

BASIS_SIZE = 8


class Preprocessing(MetaLayer):
    def __init__(self, output_dim=BASIS_SIZE, trainable=True, flav_info=None, seed=0, **kwargs):
        self.output_dim = output_dim
        self.trainable = trainable
        self.flav_info = flav_info
        self.seed = seed
        super(MetaLayer, self).__init__(**kwargs)

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.output_dim)

    def build(self, input_shape):
        self.kernel = []
        for i in range(BASIS_SIZE * 2):
            fl = int(i / 2)
            f = self.flav_info[fl]["fl"]
            if i % 2:
                name = "beta_{0}".format(f)
                kind = "largex"
            else:
                name = "alfa_{0}".format(f)
                kind = "smallx"

            limits = self.flav_info[fl][kind]
            lm = limits[0]
            lp = limits[1]
            init = MetaLayer.select_initializer("random_uniform", minval=lm, maxval=lp, seed=self.seed + i)
            const = self.constraint_MinMaxWeight(min_value=lm, max_value=lp)

            nw = self.builder_helper(
                name=name, kernel_shape=(1,), initializer=init, trainable=self.trainable, constraint=const
            )

            self.kernel.append(nw)

        super(Preprocessing, self).build(input_shape)

    def call(self, x):
        pdf_raw = self.concatenate(
            [
                x ** (1 - self.kernel[0][0]) * (1 - x) ** self.kernel[1][0],  # sigma
                x ** (1 - self.kernel[2][0]) * (1 - x) ** self.kernel[3][0],  # g
                x ** (1 - self.kernel[4][0]) * (1 - x) ** self.kernel[5][0],  # v
                x ** (1 - self.kernel[6][0]) * (1 - x) ** self.kernel[7][0],  # v3
                x ** (1 - self.kernel[8][0]) * (1 - x) ** self.kernel[9][0],  # v8
                x ** (1 - self.kernel[10][0]) * (1 - x) ** self.kernel[11][0],  # t3 = sigma
                x ** (1 - self.kernel[12][0]) * (1 - x) ** self.kernel[13][0],  # t8 = sigma
                x ** (1 - self.kernel[14][0]) * (1 - x) ** self.kernel[15][0],  # t15 c-
            ],
            axis=1,
        )
        return pdf_raw


class Rotation(MetaLayer):
    def __init__(self, output_dim=14, **kwargs):
        self.output_dim = output_dim
        super(MetaLayer, self).__init__(**kwargs, name="evolution")

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.output_dim)

    def call(self, x_raw):
        x = self.transpose(x_raw)
        pdf_raw_list = [
            0 * x[0],  # photon
            x[0],  # sigma
            x[1],  # g
            x[2],  # v
            x[3],  # v3
            x[4],  # v8
            x[2],  # v15
            x[2],  # v24
            x[2],  # v35
            x[5],  # t3
            x[6],  # t8
            x[0] - 4 * x[7],  # t15 (c-)
            x[0],  # t24
            x[0],  # t35
        ]
        pdf = self.concatenate(pdf_raw_list, target_shape=(self.output_dim, -1))
        return self.transpose(pdf)
