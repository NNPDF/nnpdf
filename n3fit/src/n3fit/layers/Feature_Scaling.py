import numpy as np
from sklearn.preprocessing import MinMaxScaler, StandardScaler

from n3fit.backends import MetaLayer


class Feature_Scaling(MetaLayer):
    """
        Applies a Normalisation of the x-grid distribution. 
    """
    
    def __init__(self, scale_features=True, **kwargs):
        self.scale_features = scale_features
        if self.scale_features:
            self.feature_range = (-1,1)
            fake_x = np.concatenate(
                (np.logspace(-6, -3, num=50, endpoint=False), np.linspace(1e-3, 1, num=50))
            )
            self.scale_ = (self.feature_range[1] - self.feature_range[0])
            self.min_ = self.feature_range[0] - fake_x.min() * self.scale_

        super().__init__(**kwargs)

    def call(self, x_raw):
        x = x_raw
        if self.scale_features:
            x *= self.scale_
            x += self.min_
        return x
