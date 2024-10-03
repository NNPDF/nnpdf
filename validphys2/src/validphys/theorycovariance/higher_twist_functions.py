import numpy as np
import pandas as pd
import scipy.linalg as la
import scipy.interpolate as scint
from validphys.convolution import central_fk_predictions

import operator

def beta_tilde_5pt(delta_h, idx):
  shifted_list = [0 for _ in range(len(delta_h))]
  shifted_list[idx] = delta_h[idx]
  return shifted_list

def ht_parametrisation(
        delta_h: list,
        nodes: list,
        x: list,
        Q2: list,
      ):
        H = scint.CubicSpline(nodes, delta_h)
        H = np.vectorize(H)

        PC = H(x) / Q2
        return PC


def null_func(size):
    """Auxiliary function used to return arrays of zeros for those datasets
    that don't require higher twist."""
    def zeros(list):
        return np.zeros(size)
    return zeros


def DIS_F2R_ht(experiment,
               pdf,
               H2p_dict,
               H2d_dict,
               x,
               q2):
  cuts = experiment.cuts
  fkspec_F2D, fkspec_F2P = experiment.fkspecs
  fk_F2D = fkspec_F2D.load_with_cuts(cuts)
  fk_F2P = fkspec_F2P.load_with_cuts(cuts)
  F2D = central_fk_predictions(fk_F2D, pdf)
  F2P = central_fk_predictions(fk_F2P, pdf)

  F2D = np.concatenate(F2D.values)
  F2P = np.concatenate(F2P.values)
  F2_ratio = operator.truediv(F2D, F2P)


  def PC_2_p(list):
      PC_p = ht_parametrisation(list, H2p_dict["nodes"], x, q2)
      result = np.array(operator.truediv(F2D, np.sum([F2P, PC_p],axis=0)) - F2_ratio)
      return result 
  
  def PC_2_d(list):
      PC_d = ht_parametrisation(list, H2d_dict["nodes"], x, q2)
      result = np.array(operator.truediv(np.sum([F2D, PC_d],axis=0), F2P) - F2_ratio)
      #result = np.array(operator.truediv(F2D, np.sum([F2D, PC_d],axis=0)) - F2_ratio) #old implementation
      return result
  
  return PC_2_p, PC_2_d


def DIS_F2_ht(H2_dict, x, q2):
  def PC_2(list):
      result = ht_parametrisation(list, H2_dict['nodes'], x, q2)
      return result
    
  return PC_2


def DIS_F2C_ht(H2p_dict,
               H2d_dict,
               x,
               q2):
  # Iron target
  Z = 23.403
  A = 49.618

  def PC_2_p(list):
    result = ht_parametrisation(list, H2p_dict['nodes'], x, q2)
    return result
  
  def PC_2_d(list):
    result = 2 * (Z - A) / A * ht_parametrisation(list, H2d_dict['nodes'], x, q2)
    return result
  
  return PC_2_p, PC_2_d


def DIS_NC_ht(H2_dict,
              HL_dict,
              x,
              q2,
              y):
  yp = 1 + np.power(1 - y, 2)
  yL = np.power(y, 2)
  N_L = - yL / yp

  def PC_2(list):
    return ht_parametrisation(list, H2_dict['nodes'], x, q2)
  
  def PC_L(list):
    return N_L * ht_parametrisation(list, HL_dict['nodes'], x, q2)
  
  return PC_2, PC_L


class DIS_SNU:
  GF2 = 1.1663787e-05 #GeV^-2
  def __init__(self,
               HT_dict,
               target_tuple, #(A,Z)
               kin_tuple,
               Mh,
               Mw,
               lepton):

    # Lead target
    A = target_tuple[0]
    Z = target_tuple[1]
    x = kin_tuple[0]
    q2 = kin_tuple[1]
    y = kin_tuple[2]
    self.nuclear_target = 2 * (Z - A) / A 
    self.Mh = Mh #0.938 GeV
    self.Mw2 = np.power(Mw, 2) # GeV^2
    self.yp = 1 + np.power(1 - y, 2) - 2 * np.power( x * y * Mh, 2) / q2
    self.yL = np.power(y, 2)
    self.ym = 1 - np.power(1 - y, 2)
    self.N = self.GF2 * Mh / ( 2 * np.pi * np.power( 1 + q2 / self.Mw2, 2) ) * self.yp
    self.H2p_dict = HT_dict['H2p']
    self.H2d_dict = HT_dict['H2d']
    self.HLp_dict = HT_dict['HLp']
    self.HLd_dict = HT_dict['HLd']
    self.H3p_dict = HT_dict['H3p']
    self.H3d_dict = HT_dict['H3d']

    self.x = x
    self.q2 = q2
    self.l = lepton
  
  def PC_2_p(self, list):
    norm = self.N
    return norm * ht_parametrisation(list, self.H2p_dict['nodes'], self.x, self.q2)
  
  def PC_2_d(self, list):
    norm = self.N * self.nuclear_target
    return norm * ht_parametrisation(list, self.H2d_dict['nodes'], self.x, self.q2)

  def PC_L_p(self, list):
    norm = - self.N * self.yL / self.yp
    return norm * ht_parametrisation(list, self.HLp_dict['nodes'], self.x, self.q2)
  
  def PC_L_d(self, list):
    norm = - self.N * self.yL / self.yp * self.nuclear_target
    return norm * ht_parametrisation(list, self.HLd_dict['nodes'], self.x, self.q2)
  
  def PC_3_p(self, list):
    norm = self.N * np.power(-1, self.l) * self.ym / self.yp * self.x
    return norm * ht_parametrisation(list, self.H3p_dict['nodes'], self.x, self.q2)
     
  def PC_3_d(self, list):
    norm = self.N * np.power(-1, self.l) * self.ym / self.yp * self.x * self.nuclear_target
    return norm * ht_parametrisation(list, self.H3d_dict['nodes'], self.x, self.q2)
  

class DIS_NUTEV(DIS_SNU):
   def __init__(self, HT_dict,
               target_tuple, #(A,Z)
               kin_tuple,
               Mh,
               Mw,
               lepton):
      super.__init__(self,HT_dict, target_tuple, kin_tuple, Mh, Mw, lepton)
      self.N = 50 * self.yp / np.power( 1 + self.q2 / self.Mw2, 2)


class DIS_HERA_CC(DIS_SNU):
   def __init__(self, HT_dict,
               kin_tuple,
               Mh,
               Mw,
               lepton):
      super.__init__(self,HT_dict, (1,1), kin_tuple, Mh, Mw, lepton)
      y = kin_tuple[2]
      self.yp = 1 + np.power(1 - y, 2)
      N = 1 / 4 * self.yp