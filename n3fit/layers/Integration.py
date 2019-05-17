from backends import MetaLayer

class xDivide(MetaLayer):
    """
    Divide the pdf by x
    By default: dimension = 8 dividing the entries corresponding to v, v3, v8
    """

    def __init__(self, output_dim = 8, div_list = [2,3,4], size = 1000, **kwargs):
        self.output_dim = output_dim
        self.div_list = div_list
        self.one = self.tensor_ones( (size, 1) )
        super(MetaLayer, self).__init__(**kwargs)

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.output_dim)

    def call(self, x):
        out_array = []
        for i in range(self.output_dim):
            if i in self.div_list:
                res = self.one/x
            else:
                res = self.one
            out_array.append(res)
        out_tensor = self.concatenate(out_array, axis = 1)
        return out_tensor

class xMultiply(MetaLayer):
    """
    Multiply the pdf by x
    By default: dimension = 8 multiply all entries but v, v3, v8
    """

    def __init__(self, output_dim = 8, div_list = [2,3,4], size = 1000, **kwargs):
        self.output_dim = output_dim
        self.div_list = div_list
        self.one = self.tensor_ones( (size, 1) )
        super(MetaLayer, self).__init__(**kwargs)

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.output_dim)

    def call(self, x):
        out_array = []
        for i in range(self.output_dim):
            if i in self.div_list:
                res = self.one
            else:
                res = self.one*x
            out_array.append(res)
        out_tensor = self.concatenate(out_array, axis = 1)
        return out_tensor

class xIntegrator(MetaLayer):
    """
    This layer performs a sum of the input layer/tensor on the first axis
    """
    def __init__(self, grid_weights, output_dim = 8, **kwargs):
        grid_weights_tensor = self.np_to_tensor(grid_weights)
        # Open up the grid weights
        self.grid_weights = self.many_replication(grid_weights_tensor, replications = 8, axis = 1)
        self.output_dim = output_dim
        super(MetaLayer, self).__init__(**kwargs)

    def compute_output_shape(self, input_shape):
        return (self.output_dim, )

    def call(self, x):
        xx = x*self.grid_weights
        return self.sum(xx, axis = 0)


class MSR_Normalization(MetaLayer):
    """
    Applies the normalisation so that the PDF output fullfills the sum rules
    """
    def __init__(self, output_dim = 14, **kwargs):
        self.output_dim = output_dim
        self.one = self.tensor_ones( (1, 1) )
        self.three = 3*self.tensor_ones( (1, 1) )
        super(MSR_Normalization, self).__init__(**kwargs, name = "normalizer")

    def compute_output_shape(self, input_shape):
        return (self.output_dim, )


    def call(self, integrated):
        x = integrated
        pdf_sr = self.concatenate([
            self.one, # photon
            self.one, # sigma
            (self.one - x[0])/x[1], # g 
            self.three/x[2], # v
            self.one/x[3], # v3
            self.three/x[4], # v8
            self.three/x[2], # v15
            self.three/x[2], # v24
            self.three/x[2], # v35
            self.one, # t3
            self.one, # t8
            self.one, # t15 (c-)
            self.one, # t24
            self.one, # t35
        ])
        return pdf_sr
