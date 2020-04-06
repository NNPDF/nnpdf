import numpy as np
from sklearn.preprocessing import MinMaxScaler, StandardScaler

from n3fit.backends import MetaLayer


class Feature_Scaling(MetaLayer):
    """
        Applies a Normalisation of the x-grid distribution. 
    """

    def __init__(self, scale_features=True, **kwargs):
        self.scale_features = scale_features
        self.first_run = True
        self.scaler = None
        super().__init__(**kwargs)

    def scale_features_func(self, x_raw):
        fake_x = np.concatenate(
            (np.logspace(-6, -3, num=50, endpoint=False), np.linspace(1e-3, 1, num=50))
        )
        self.scaler = MinMaxScaler(feature_range=(0.5, 1), copy=False)
        self.scaler.fit(fake_x.reshape(-1,1))

    def call(self, x_raw):
        x = x_raw
        if self.scale_features:
            if self.first_run:
                self.first_run = False
                self.scale_features_func(x_raw=x)
            self.scaler.transform(x.tensor_content.reshape(-1, 1))
        return x
