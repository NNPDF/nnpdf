# -*- coding: utf-8 -*-
"""
kintransforms.py

This modules defines classes that transform the kinematics as defined in the
CommonData files into some function of these kinematics that is more convenient
for representation.The kinematic transforms should also define an 'xq2map',
that maps each kinematic point into zero or more points in (x, Q²), as a
function of the **new** kinematics.

The expected interface of the classes is:

    .. code-block:: python

        class mytransform:
            def __call__(self, k1:np.array,k2:np.array,k3:np.array) -> (np.array, np.array, np.array):
                #Transform kinematics
                ...
                return trasformed_k1, transformed_k2, transformed_k3

            def new_labels(self, old_label1:str, old_label2:str, old_label3:str) -> (str, str, str):
                #Transform labels
                return transformed_label1, transformed_label2, transformed_label3

            #Using as input the result of __call__ as well as any new labels
            def xq2map(self, k1:np.array,k2:np.array,k3:np.array,**extra_labels) -> (np.array, np.array):
                #calculate (x,Q²)
                return x, Q2


"""

#TODO: fix the issue with Zmass and top mass - make them (globa?) constants
ZMASS=91.1876
TMASS=173.3

import abc

import numpy as np

class Kintransform(metaclass=abc.ABCMeta):
    @classmethod
    def __subclasshook__(cls, other):
        return hasattr(other, 'xq2map') and hasattr(other, '__call__') and hasattr(other, 'new_labels')

#Common utilities on top of which we build the transforms

class SqrtScaleMixin:
    def __call__(self, k1, k2, k3):
        return  k1, np.sqrt(k2), k3

    qlabel = NotImplemented

    def new_labels(self, s1, s2, s3):
        return s1, self.qlabel, s3

class DISXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        """in DIS-like experiment k1 is x, k2 is Q"""
        return k1, k2*k2

class DYXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        """in DY-like experiments k1 is (pseudo)-rapidity and k2 is Q for each point in the experiment there are two points in the xQ2 map"""
        ratio = k2/k3
        x1 = ratio*np.exp(k1)
        x2 = ratio*np.exp(-k1)
        q2 = k2*k2
        x = np.concatenate(( x1,x2 ))
        return np.clip(x,None,1, out=x), np.concatenate(( q2,q2 ))

class JETXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        """in DY-like experiments k1 is (pseudo)-rapidity and k2 is pT"""
        ratio = k2/k3
        x = ratio*(np.exp(k1)+np.exp(-k1))
        q2 = k2*k2
        return np.clip(x,None,1, out=x), q2

class EWPTXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        """in ZPt-like Experiments k1 is the pt, k2 is Q"""
        zmass2 = ZMASS*ZMASS
        Q = (np.sqrt(zmass2+k1*k1)+k1)
        effQ = np.sqrt(zmass2+k1*k1)
        return Q/k3, effQ*effQ

class DYMXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        """in DYM-like experiments the k1 is the mass, k2 is the mass"""
        return k2/k3, k2*k2

class HQPTXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        """in HQPt-like Experiments k1 is the pt, k2 is Q"""
        QMASS2 = TMASS*TMASS
        Q = (np.sqrt(QMASS2+k1*k1)+k1)
        return Q/k3, Q*Q

class HQQPTXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        """in ZPt-like Experiments k1 is the pt, k2 is Q"""
        QQMASS2 = (2*TMASS)*(2*TMASS)
        Q = (np.sqrt(QQMASS2+k1*k1)+k1)
        return Q/k3, Q*Q

#The transforms themselves
class identity:
    def __call__(self, k1, k2, k3):
        return k1, k2, k3

    def new_labels(self, *old_labels):
        return old_labels

class dyp_sqrt_scale(SqrtScaleMixin, DYXQ2MapMixin):
    qlabel = '$M (GeV)$'


class jet_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin):
    def new_labels(self, *old_labels):
        return ('$|y|$', '$p_T$ (GeV)', r'$\sqrt{s} (GeV)$')

class dis_sqrt_scale(DISXQ2MapMixin):
    def __call__(self, k1, k2, k3):
        ecm = np.sqrt(k2/(k1*k3))
        return k1, np.sqrt(k2), np.ceil(ecm)

    def new_labels(self, *old_labels):
        return ('$x$', '$Q$ (GeV)', r'$\sqrt{s} (GeV)$')

class ewj_jpt_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin):
    qlabel = '$M (GeV)$'

class ewj_jrap_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin):
    qlabel = '$M (GeV)$'

class ewj_mll_sqrt_scale(SqrtScaleMixin,DYMXQ2MapMixin):
    qlabel = '$M_{ll} (GeV)$'

class ewj_pt_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin):
    qlabel = '$M (GeV)$'

class ewj_ptrap_sqrt_scale(SqrtScaleMixin,EWPTXQ2MapMixin):
    qlabel = r'$p_T (GeV)$'

class ewj_rap_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin):
    qlabel = '$M (GeV)$'

class ewk_mll_sqrt_scale(SqrtScaleMixin,DYMXQ2MapMixin):
    qlabel = '$M_{ll} (GeV)$'

class ewk_pt_sqrt_scale(SqrtScaleMixin,EWPTXQ2MapMixin):
    qlabel = '$M (GeV)$'

class ewk_ptrap_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin):
    qlabel = r'$p_T (GeV)$'

class ewk_rap_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin):
    qlabel = '$M (GeV)$'

class hig_rap_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin):
    qlabel = '$M_H (GeV)$'

class hqp_mqq_sqrt_scale(SqrtScaleMixin,DYMXQ2MapMixin):
    qlabel = r'$\mu (GeV)$'

class hqp_ptq_sqrt_scale(SqrtScaleMixin,HQPTXQ2MapMixin):
    qlabel = r'$\mu (GeV)$'

class hqp_ptqq_sqrt_scale(SqrtScaleMixin,HQQPTXQ2MapMixin):
    qlabel = r'$\mu (GeV)$'

class hqp_yq_sqrt_scale(SqrtScaleMixin,JETXQ2MapMixin):
    qlabel = r'$\mu (GeV)$'

class hqp_yqq_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin):
    qlabel = r'$\mu (GeV)$'

class inc_sqrt_scale(SqrtScaleMixin,DYMXQ2MapMixin):
    qlabel = r'$\mu (GeV)$'

class pht_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin):
    qlabel = r'$E_{T,\gamma} (GeV)$'

class sia_sqrt_scale(SqrtScaleMixin,DISXQ2MapMixin):
    qlabel = '$Q (GeV)$'


class nmc_process(DISXQ2MapMixin):
    def __call__(self, k1,k2,k3):
        xBins = [0.0045, 0.008, 0.0125, 0.0175,
                 0.025, 0.035, 0.05, 0.07, 0.09, 0.11,
                 0.14, 0.18, 0.225, 0.275, 0.35, 0.5]
        for x in np.nditer(k1, op_flags=['readwrite']):
            x[...] = min(xBins, key=lambda y:abs(x-y))
        ecm = np.sqrt(k2/(k1*k3))
        return k1, np.sqrt(k2), np.ceil(ecm)

    def new_labels(self, *old_labels):
        return ('$x$', '$Q$ (GeV)', r'$\sqrt{s} (GeV)$')