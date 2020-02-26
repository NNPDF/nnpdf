from n3fit.backends import MetaLayer


class NewPositivity(MetaLayer):
    """
        Just a Template for now

        This layer is just for testing so it is not necessariliy efficient
        the idea is that it takes as input
            (pdf(x), pdf(z))
        where pdf is the positivity-scheme pdf and x is the xgrid for the target fktable
        and returns a pdf(x) in the MSbar scheme

        The x is given at instantiation time in the variable x_in.
        For this first test we take the x and the z to be the same grid,
        so life is easier
    """

    def __init__(self, x_in = None, **kwargs):
        self.x_in = x_in
        super().__init__(**kwargs)

    def call(self, x):
        pdf_in_x = x[0] # tensor with indices (i, f) where i is the index of the grid in x and f the flavour index
        pdf_in_z = x[1]
        # Very non-trivial operation that changes completely the PDF happens here
        result = pdf_in_z + pdf_in_x
        return result
