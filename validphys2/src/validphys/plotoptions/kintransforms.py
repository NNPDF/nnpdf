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


The kinematic labels are:

    .. code-block:: python

        {'DIS': ('$x$', '$Q^2 (GeV^2)$', '$y$'),
         'DYP': ('$y$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'EWJ_JPT': ('$p_T (GeV)$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'EWJ_JRAP': ('$\\eta/y$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'EWJ_MLL': ('$M_{ll} (GeV)$', '$M_{ll}^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'EWJ_PT': ('$p_T (GeV)$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'EWJ_PTRAP': ('$\\eta/y$', '$p_T^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'EWJ_RAP': ('$\\eta/y$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'EWK_MLL': ('$M_{ll} (GeV)$', '$M_{ll}^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'EWK_PT': ('$p_T$ (GeV)', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'EWK_PTRAP': ('$\\eta/y$', '$p_T^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'EWK_RAP': ('$\\eta/y$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'HIG_RAP': ('$y$', '$M_H^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'HQP_MQQ': ('$M^{QQ} (GeV)$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'HQP_PTQ': ('$p_T^Q (GeV)$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'HQP_PTQQ': ('$p_T^{QQ} (GeV)$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'HQP_YQ': ('$y^Q$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'HQP_YQQ': ('$y^{QQ}$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'INC': ('$0$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'JET': ('$\\eta$', '$p_T^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'DIJET': ('$\\eta$', '$m_{12} (GeV)$', '$\\sqrt{s} (GeV)$'),
         'PHT': ('$\\eta_\\gamma$', '$E_{T,\\gamma}^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
         'SIA': ('$z$', '$Q^2 (GeV^2)$', '$y$')}



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
        """in DY-like experiments k1 is (pseudo)-rapidity and k2 is Q for
        each point in the experiment there are two points in the xQ2 map"""
        ratio = k2/k3
        x1 = ratio*np.exp(k1)
        x2 = ratio*np.exp(-k1)
        q2 = k2*k2
        x = np.concatenate(( x1,x2 ))
        return np.clip(x,a_min=None,a_max=1, out=x), np.concatenate(( q2,q2 ))

class JETXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        """
        k1 is (pseudo)-rapidity and k2 is pT
        plotting both x1 and x2
        """
        ratio = k2/k3
        x1 = 2 * ratio * np.exp(k1)
        x2 = 2 * ratio * np.exp(-k1)
        q2 = k2*k2
        x = np.concatenate(( x1,x2 ))
        return np.clip(x,a_min=None,a_max=1, out=x), np.concatenate(( q2,q2 ))

class DIJETXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        """
        k1 is max(|y1|,|y2|) and k2 is m12
        plotting both x1 and x2
        """
        ratio = k2/k3
        x1 = ratio * np.exp(k1)
        x2 = ratio * np.exp(-k1)
        q2 = k2*k2
        x = np.concatenate(( x1,x2 ))
        return np.clip(x,a_min=None,a_max=1, out=x), np.concatenate(( q2,q2 ))

class DIJETATLASXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        """
        k1 is rapidity difference and k2 is m12
        plotting both x1 and x2
        """
        ratio = k2/k3
        #x1 = ratio
        #x2 = np.full_like(x1, 1.0)
        x1 = ratio * np.exp(k1)
        x2 = ratio * np.exp(-k1)
        q2 = k2*k2
        x = np.concatenate(( x1,x2 ))
        return np.clip(x,a_min=None,a_max=1, out=x), np.concatenate(( q2,q2 ))
    
class DIJET3DXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        """
        k1 is the rapidity difference, k2 is pTavg, k3 is boost rapidity
        TODO: NB!! HARDCODING sqrt(s) = 8 TeV, the c.m. energy of 1705.02628.pdf
        plotting both x1 and x2
        """
        sqrts = 8000
        ratio = k2/sqrts
        prefactor = ratio * (np.exp(k1) + np.exp(-k1))
        x1 = prefactor * np.exp(k3)
        x2 = prefactor * np.exp(-k3)
        q2 = k2*k2
        x = np.concatenate(( x1,x2 ))
        #print(k1[55],k2[55],k3[55],x[55],np.argwhere(x>1))
        return np.clip(x,a_min=None,a_max=1, out=x), np.concatenate(( q2,q2 ))



    
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


class dyp_sqrt_scale(SqrtScaleMixin, DYXQ2MapMixin):
    qlabel = '$M (GeV)$'


class jet_sqrt_scale(SqrtScaleMixin,JETXQ2MapMixin):
    def new_labels(self, *old_labels):
        return ('$|y|$', '$p_T$ (GeV)', r'$\sqrt{s} (GeV)$')

class dijet_sqrt_scale(SqrtScaleMixin,DIJETXQ2MapMixin):
    def new_labels(self, *old_labels):
        return ('$|y|$', '$m_{12}$ (GeV)', r'$\sqrt{s} (GeV)$')

class dijet_sqrt_scale_ATLAS(SqrtScaleMixin,DIJETATLASXQ2MapMixin):
    def __call__(self, k1, k2, k3):
        return k1, k2, k3

    def new_labels(self, *old_labels):
        return ('$|y^*|$', '$m_{12}$ (GeV)', r'$\sqrt{s} (GeV)$')

class dijet_CMS_3D(SqrtScaleMixin,DIJET3DXQ2MapMixin):
    def new_labels(self, *old_labels):
        return ('$|y^*|$', '$p_{T,avg}$ (GeV)', r'$|y_b|$')

class dijet_CMS_5TEV(SqrtScaleMixin,DIJET3DXQ2MapMixin):
    def new_labels(self, *old_labels):
        return ('$\eta_{dijet}$', '$p_{T,avg}$ (GeV)', r'$\sqrt{s} (GeV)$')    

class dis_sqrt_scale(DISXQ2MapMixin):
    def __call__(self, k1, k2, k3):
        ecm = np.sqrt(k2/(k1*k3))
        return k1, np.sqrt(k2), np.ceil(ecm)

    def new_labels(self, *old_labels):
        return ('$x$', '$Q$ (GeV)', r'$\sqrt{s} (GeV)$')

class ewj_jpt_sqrt_scale(SqrtScaleMixin,EWPTXQ2MapMixin):  #okay but it does not exist
    qlabel = '$M (GeV)$'

class ewj_jrap_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin): #EWJ_JRAP->DY  ----> okay but it does not exist
    qlabel = '$M (GeV)$'

class ewj_mll_sqrt_scale(SqrtScaleMixin,DYMXQ2MapMixin): #EWJ_MLL->DYm  ----> okay but it does not exist
    qlabel = '$M_{ll} (GeV)$'

class ewj_pt_sqrt_scale(SqrtScaleMixin,EWPTXQ2MapMixin): #EWJ_PT->DY ----> Zpt, okay but it does not exist
    qlabel = '$M (GeV)$'

class ewj_ptrap_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin): # EWJ_PTRAP -> DY okay, but it does not exist
    qlabel = r'$p_T (GeV)$'

class ewj_rap_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin):  # EWJ_RAP -> DY okay (can we get rid of it also in commondata?)
    qlabel = '$M (GeV)$'

class ewk_mll_sqrt_scale(SqrtScaleMixin,DYMXQ2MapMixin): # EWK_MLL -> DYM okay
    qlabel = '$M_{ll} (GeV)$'

class ewk_pt_sqrt_scale(SqrtScaleMixin,EWPTXQ2MapMixin): # EWK_PT -> Zpt okay
    qlabel = '$M (GeV)$'

class ewk_ptrap_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin): # EWK_PT -> DY okay
    qlabel = r'$p_T (GeV)$'

class ewk_rap_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin): # EWK_RAP -> DY okay
    qlabel = '$M (GeV)$'

class hig_rap_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin): #okay, but it does not exist
    qlabel = '$M_H (GeV)$'

class hqp_mqq_sqrt_scale(SqrtScaleMixin,DYMXQ2MapMixin): # HQP_MQQ -> DYM okay
    qlabel = r'$\mu (GeV)$'

class hqp_ptq_sqrt_scale(SqrtScaleMixin,HQPTXQ2MapMixin): # HQP_PTQ -> HQPT okay
    qlabel = r'$\mu (GeV)$'

class hqp_ptqq_sqrt_scale(SqrtScaleMixin,HQQPTXQ2MapMixin): # HQP_PTQQ -> HQQPT okay
    qlabel = r'$\mu (GeV)$'

class hqp_yq_sqrt_scale(SqrtScaleMixin,JETXQ2MapMixin): # HQP_YQ->JETXQ2 okay
    qlabel = r'$\mu (GeV)$'

class hqp_yqq_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin): #HQP_YQQ->DYXQ2 okay
    qlabel = r'$\mu (GeV)$'

class inc_sqrt_scale(SqrtScaleMixin,DYMXQ2MapMixin): # INC -> DYM okay
    qlabel = r'$\mu (GeV)$'

class pht_sqrt_scale(SqrtScaleMixin,DYXQ2MapMixin): #okay but not in commondata
    qlabel = r'$E_{T,\gamma} (GeV)$'

class sia_sqrt_scale(SqrtScaleMixin,DISXQ2MapMixin): #okay but not in commondata
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

class ewk_pseudorapity_sqrt_scale(ewk_rap_sqrt_scale):
    def new_labels(self, *old_labels):
        superlabels = super().new_labels(*old_labels)
        return (r'$\eta$', *superlabels[1:])
