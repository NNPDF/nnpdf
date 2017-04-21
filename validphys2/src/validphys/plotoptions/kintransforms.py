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
import abc

from numpy import sqrt, ceil, nditer

class Kintransform(metaclass=abc.ABCMeta):
    @classmethod
    def __subclasshook__(cls, other):
        return hasattr(other, '__call__') and hasattr(other, 'new_labels')

#TODO: Do we ever allow a transform without this?
class KintransformWithXQ2Map(Kintransform):
    @classmethod
    def __subclasshook__(cls, other):
        return hasattr(other, 'xq2map') and issubclass(other, Kintransform)

#Common utilities on top of which we build the transforms

class SqrtScaleMixin:
    def __call__(self, k1, k2, k3):
        return  k1, sqrt(k2), k3

class DISXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        #in DIS-like experiment k1 is x, k2 is Q 
        return k1, k2*k2

class DYXQ2MapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        #in DY-like experiments k1 is (pseudo)-rapidity and k2 is Q
        #for each point in the experiment there are 
        #two points in the xQ2 map
        x1 = k2/k3*Exp(k1)
        x2 = k2/k3*Exp(-k1)
        return np.concatenate(( x1,x2 )), np.concatenate(( k2*k2,k2*k2 ))

class ZPTMapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        #in ZPt-like Experiments k1 is the pt, k2 is Q
        zmass = 91.1876
        Q = (np.sqrt(zmass^2+k1*k1)+k1)
        return Q/k3, Q*Q

class JetPTMapMixin:
    def xq2map(self, k1, k2, k3, **extra_labels):
        #in JetPt-like Experiments k1 is the pt, k2 is Q
        return k1/k3, k1*k1

#The transforms themselves
class identity:
    def __call__(self, k1, k2, k3):
        return k1, k2, k3

    def new_labels(self, *old_labels):
        return old_labels

class dyp_sqrt_scale(SqrtScaleMixin):
    def new_labels(self, *old_labels):
        return ('$y$', '$M$ (GeV)', r'$\sqrt{s} (GeV)$')

class jet_sqrt_scale(SqrtScaleMixin):
    def new_labels(self, *old_labels):
        return ('$|y|$', '$p_T$ (GeV)', r'$\sqrt{s} (GeV)$')

class dis_sqrt_scale(DISXQ2MapMixin):
    def __call__(self, k1, k2, k3):
        ecm = sqrt(k2/(k1*k3))
        return k1, sqrt(k2), ceil(ecm)

    def new_labels(self, *old_labels):
        return ('$x$', '$Q$ (GeV)', r'$\sqrt{s} (GeV)$')

class nmc_process(DISXQ2MapMixin):
    def __call__(self, k1,k2,k3):
        xBins = [0.0045, 0.008, 0.0125, 0.0175,
                 0.025, 0.035, 0.05, 0.07, 0.09, 0.11,
                 0.14, 0.18, 0.225, 0.275, 0.35, 0.5]
        for x in nditer(k1, op_flags=['readwrite']):
            x[...] = min(xBins, key=lambda y:abs(x-y))
        ecm = sqrt(k2/(k1*k3))
        return k1, sqrt(k2), ceil(ecm)

    def new_labels(self, *old_labels):
        return old_labels