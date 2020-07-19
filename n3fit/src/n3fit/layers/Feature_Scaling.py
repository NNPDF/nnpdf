import numpy as np
from functools import reduce

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

fk_global_datasets = [
    "NMCPD_D",
    # 'NMCPD_P',
    "NMC",
    "SLACP",
    "SLACD",
    "BCDMSP",
    "BCDMSD",
    "CHORUSNU",
    "CHORUSNB",
    "NTVNUDMN",
    "NTVNBDMN",
    "HERACOMBNCEM",
    "HERACOMBNCEP460",
    "HERACOMBNCEP575",
    "HERACOMBNCEP820",
    "HERACOMBNCEP920",
    "HERACOMBCCEM",
    "HERACOMBCCEP",
    "HERAF2CHARM",
    "H1HERAF2B",
    "ZEUSHERAF2B",
    # "DYE886R_D",
    "DYE886R_P",
    "DYE886P",
    "DYE605",
    "CDFZRAP",
    "CDFR2KT",
    "D0ZRAP",
    "D0WEASY_WM",
    # 'D0WEASY_WP',
    "D0WMASY_WM",
    # 'D0WMASY_WP',
    "ATLASWZRAP36PB",
    "ATLASZHIGHMASS49FB",
    "ATLASLOMASSDY11EXT",
    "ATLASWZRAP11",
    "ATLASR04JETS36PB",
    "ATLASR04JETS2P76TEV",
    "ATLAS1JET11",
    "ATLASZPT8TEVMDIST",
    "ATLASZPT8TEVYDIST",
    "ATLASTTBARTOT",
    "ATLASTOPDIFF8TEVTRAP_NUM",
    # 'ATLASTOPDIFF8TEVTRAP_DEN',
    "CMSWEASY840PB_WM",
    # 'CMSWEASY840PB_WP',
    "CMSWCHARM_WP",
    # 'CMSWCHARM_WM',
    "CMSWMU8TEV",
    "CMSJETS11",
    "CMS1JET276TEV",
    "CMSZDIFF12",
    "CMSTTBARTOT",
    "CMSTOPDIFF8TEVTTRAP_NUM",
    # 'CMSTOPDIFF8TEVTTRAP_DEN',
    "LHCBZ940PB",
    "LHCBZEE2FB",
    "LHCBWZMU7TEV",
    "LHCBWZMU8TEV",
]

l = Loader() 

class Feature_Scaling(MetaLayer):
    """
        Applies a Normalisation of the x-grid distribution.
    """

    def __init__(
        self, **kwargs
    ):
        fk_xgrids = np.array([])
        for fk_dataset in fk_dis_datasets:
            print(f"loading {fk_dataset}")
            fk = l.check_fktable(setname=fk_dataset, theoryID=53, cfac=[])
            res = load_fktable(fk)
            fk_xgrids = np.concatenate([fk_xgrids, res.xgrid])
        fk_xgrids = np.sort(fk_xgrids)

        xgrid_size = fk_xgrids.size
        new_xgrid = np.linspace(start=-1, stop=1, num=xgrid_size)
        unique, counts = np.unique(fk_xgrids, return_counts=True)

        self.map_from = unique
        self.map_to = []
        for cumsum_ in np.cumsum(counts):
            self.map_to.append(new_xgrid[cumsum_ - 1])
        self.map_to = np.array(self.map_to)

        # map_to_temp = []
        # map_from_temp = []
        # for num in range(len(self.map_to)):
        #     if num%50==0:
        #         map_to_temp.append(self.map_to[num])
        #         map_from_temp.append(self.map_from[num])
        # self.map_to = np.array(map_to_temp)
        # self.map_from = np.array(map_from_temp)

        super().__init__(**kwargs)

    def call(self, xgrid_in):

        smallx_cond = tf.math.less(xgrid_in, self.map_to.min())
        largex_cond = tf.math.greater(xgrid_in, self.map_to.max())

        inter_cond = []
        y_inter_coefs = []
        for num in range(len(self.map_from)-1):
            x1 = self.map_from[num]
            x2 = self.map_from[num+1]
            y1 = self.map_to[num]
            y2 = self.map_to[num+1]

            y_inter_coefs.append( [x1,x2,y1,y2] )

            geq_cond = tf.math.greater_equal(xgrid_in, x1)
            less_cond = tf.math.less(xgrid_in, x2)

            inter_cond.append( tf.math.logical_and(geq_cond, less_cond) )

        def y_inter(x, coefs):
            return coefs[2] + (x - coefs[0]) * (coefs[2]-coefs[3])/(coefs[0]-coefs[1])

        res = xgrid_in
        for num in range( len(self.map_from)-1 ):        
            res = tf.where(inter_cond[num], y_inter(res, y_inter_coefs[num]), res)

        res = tf.where(smallx_cond, y_inter(res, y_inter_coefs[0]), res)
        res = tf.where(largex_cond, y_inter(res, y_inter_coefs[-1]), res)

        return res

    # def call(self, xgrid_in):

    #     smallx_cond = tf.math.less(xgrid_in, self.map_to.min())
    #     largex_cond = tf.math.greater_equal(xgrid_in, self.map_to.max())

    #     inter_cond = []
    #     y_inter_coefs = []
    #     for num in range(len(self.map_from)-1):
    #         x1 = tf.constant(self.map_from[num], dtype=xgrid_in.dtype)
    #         x2 = tf.constant(self.map_from[num+1], dtype=xgrid_in.dtype)
    #         y1 = tf.constant(self.map_to[num], dtype=xgrid_in.dtype)
    #         y2 = tf.constant(self.map_to[num+1], dtype=xgrid_in.dtype)

    #         y_inter_coefs.append( [x1,x2,y1,y2] )

    #         geq_cond = tf.math.greater_equal(xgrid_in, x1)
    #         less_cond = tf.math.less(xgrid_in, x2)

    #         inter_cond.append( tf.math.logical_and(geq_cond, less_cond) )

    #     @tf.function
    #     def y_inter(x, x1, x2, x3, x4):
    #         return y1 + (x - x1) * (y1-y2)/(x1-x2)

    #     @tf.function
    #     def where_func(cond, if_true, if_false):
    #         return tf.where(cond, if_true, if_false)
        
    #     res = where_func(inter_cond[0], y_inter(xgrid_in, y_inter_coefs[0][0], y_inter_coefs[0][1], y_inter_coefs[0][2], y_inter_coefs[0][3]), xgrid_in)
    #     @tf.function
    #     def aaa(res, inter_cond, y_inter_coefs):
    #         for num in range( 1, len(self.map_from)-1 ):        
    #             res = where_func(inter_cond[num], y_inter(res, y_inter_coefs[num][0], y_inter_coefs[num][1], y_inter_coefs[num][2], y_inter_coefs[num][3]), res)
    #         return res
    #     res = aaa(res, inter_cond, y_inter_coefs)

    #     res = where_func(smallx_cond, y_inter(res, y_inter_coefs[0][0], y_inter_coefs[0][1], y_inter_coefs[0][2], y_inter_coefs[0][3]), res)
    #     res = where_func(largex_cond, y_inter(res, y_inter_coefs[-1][0], y_inter_coefs[-1][1], y_inter_coefs[-1][2], y_inter_coefs[-1][3]), res)

    #     return res
