import numpy as np

from n3fit.backends import MetaLayer
from tensorflow.keras import backend as K
import tensorflow as tf

from validphys.fkparser import load_fktable
from validphys.loader import FallbackLoader as Loader

fk_dis_datasets=[    
    'NMCPD_D',
    # 'NMCPD_P',
    'NMC',
    'SLACP',
    'SLACD',
    'BCDMSP',
    'BCDMSD',
    'CHORUSNU',
    'CHORUSNB',
    'NTVNUDMN',
    'NTVNBDMN',
    'HERACOMBNCEM',
    'HERACOMBNCEP460',
    'HERACOMBNCEP575',
    'HERACOMBNCEP820',
    'HERACOMBNCEP920',
    'HERACOMBCCEM',
    'HERACOMBCCEP',
    'HERAF2CHARM',
    'H1HERAF2B',
    'ZEUSHERAF2B',
    ]

class Feature_Scaling(MetaLayer):
    """
        Applies a Normalisation of the x-grid distribution.
    """

    def __init__(
        self, scaler=False, **kwargs
    ):

        scaler_list = ["MinMaxScaler", "GaussianScaler"]

        # name scalar type
        self.scaler = scaler
        if self.scaler != False:
            if self.scaler not in scaler_list:
                raise ValueError(f"Feature Scaling does not support '{self.scaler}'.")

        def load_fk_tables(datasets):
            l = Loader() 
            fk_xgrids = np.array([])
            for fk_dataset in datasets:
                print(f'loading {fk_dataset}')
                fk = l.check_fktable(setname=fk_dataset, theoryID=52, cfac=[])
                res = load_fktable(fk)
                fk_xgrids = np.concatenate([fk_xgrids, res.xgrid]) 
            return fk_xgrids

        if self.scaler == 'MinMaxScaler':
            scaled_xgrids = load_fk_tables(fk_dis_datasets)
            self.max_ /= scaled_xgrids.max()

        
        elif self.scaler == 'GaussianScaler':
            scaled_xgrids = load_fk_tables(fk_dis_datasets)
            min_ = scaled_xgrids.min()
            self.max_ = scaled_xgrids.max()
            self.mean_ = scaled_xgrids.mean()
            self.scale_ = (1- -1)/(self.max_ - min_)

        super().__init__(**kwargs)

    def call(self, x_raw):
        scaled_xgrids = x_raw

        def log10(x):
            numerator = K.log(x)
            denominator = K.log( tf.constant(10, dtype=numerator.dtype))
            return numerator/denominator

        if self.scaler == "MinMaxScaler":
            scaled_xgrids = scaled_xgrids + scaled_xgrids**0.4 + scaled_xgrids**0.3 + 0.5*scaled_xgrids**0.2
            scaled_xgrids = self.max_
            scaled_xgrids *= 2
            scaled_xgrids -= 1

        elif self.scaler == "GaussianScaler":
            scaled_xgrids = scaled_xgrids + scaled_xgrids**0.4 + scaled_xgrids**0.3 + 0.5*scaled_xgrids**0.2
            
            scaled_xgrids /= self.max_
            scaled_xgrids *= 2
            scaled_xgrids -= 1
            
            scaled_xgrids = K.exp((scaled_xgrids-self.mean_)**2 ) * scaled_xgrids
                
        return scaled_xgrids
