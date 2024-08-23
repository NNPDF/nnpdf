"""
    DIS layer

    This layer produces a DIS observable, which can consists of one or more fktables.
    The rationale behind this layer is to keep all required operation in one single place
    such that is easier to optimize or modify.

    Comment on the branching based on the number of replicas:
        This is purely for performance, masking the PDF is more efficient than padding the fk table
        for one replica, and so is tensordot over einsum.

        Some timings done on snellius using tensorflow 2.15.0 and varying these 2 factors:
            | CPU\GPU | einsum | tensordot |
            | -- | -- | -- |
            | mask pdf | -  | 92 \ 65 |
            |mask fk | 330 \ 53 \  | 177 \ 53 |

            These timings are all for one replica.

            Crucially, `einsum` is a requirement of the multireplica case, while `tensordot` gives a benefit of a factor of 2x for the single replica case.
            Since this branching is required anyhow,
             by masking the PDF for 1 replica instead of padding the fktable we get an extra factor of x2
"""

import numpy as np
from scipy import interpolate as scint

from n3fit.backends import operations as op

from validphys.theorycovariance.construction import compute_normalisation_by_experiment

from .observable import Observable


class DIS(Observable):
    """
    The DIS class receives a list of active flavours and a fktable
    and prepares a layer that performs the convolution of said fktable with
    the incoming pdf.

    The fktable is expected to be rank 3 (ndata, xgrid, flavours)
    while the input pdf is rank 4 of shape (batch_size, replicas, xgrid, flavours)
    """

    def __init__(self, fktable_data, fktable_arr, dataset_name, boundary_condition=None, operation_name="NULL", nfl=14, n_replicas=1, exp_kinematics=None, **kwargs):
        super().__init__(fktable_data, fktable_arr, dataset_name, boundary_condition, operation_name, nfl, n_replicas, **kwargs)

        self.power_corrections = None
        if exp_kinematics is not None:
          self.exp_kinematics = exp_kinematics
          self.power_corrections = self.compute_abmp_parametrisation()

    def compute_abmp_parametrisation(self):
        """
        This function is very similar to `compute_ht_parametrisation` in
        validphys.theorycovariance.construction.py. However, the latter
        accounts for shifts in the 5pt prescription. As of now, this function
        is meant to work only for DIS NC data, using the ABMP16 result.
        """
        x_knots = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
        y_h2 = [0.023, -0.032, -0.005, 0.025, 0.051, 0.003, 0.0]
        y_ht = [-0.319, -0.134, -0.052, 0.071, 0.030, 0.003, 0.0]
        h2_sigma = [0.019, 0.013, 0.009, 0.006, 0.005, 0.004]
        ht_sigma = [0.126, 0.040, 0.030, 0.025, 0.012, 0.007]
        H_2 = scint.CubicSpline(x_knots, y_h2)
        H_T = scint.CubicSpline(x_knots, y_ht)

        # Reconstruct HL from HT and H2
        def H_L(x):
            return (H_2(x) - np.power(x, 0.05) * H_T(x))

        H_2 = np.vectorize(H_2)
        H_L = np.vectorize(H_L)

        x = self.exp_kinematics['kin1']
        y = self.exp_kinematics['kin3']
        Q2 = self.exp_kinematics['kin2']
        N2, NL = compute_normalisation_by_experiment(self.dataname, x, y, Q2)

        PC_2 = N2 * H_2(x) / Q2
        PC_L = NL * H_L(x) / Q2
        power_correction = PC_2 + PC_L
        power_correction = power_correction.to_numpy()

        return power_correction


    def gen_mask(self, basis):
        """
        Receives a list of active flavours and generates a boolean mask tensor

        Parameters
        ----------
            basis: list(int)
                list of active flavours

        Returns
        -------
            mask: tensor
                rank 1 tensor (flavours)
        """
        if basis is None:
            self.basis = np.ones(self.nfl, dtype=bool)
        else:
            basis_mask = np.zeros(self.nfl, dtype=bool)
            for i in basis:
                basis_mask[i] = True
        return op.numpy_to_tensor(basis_mask, dtype=bool)

    def pad_fk(self, fk, mask):
        """
        Combine an fk table and a mask into an fk table padded with zeroes for the inactive
        flavours, to be contracted with the full PDF.

        Parameters
        ----------
            fk: tensor
                FK table of shape (ndata, active_flavours, x)
            mask: tensor
                mask of shape (flavours, active_flavours)

        Returns
        -------
            padded_fk: tensor
                masked fk table of shape ndata, x, flavours)
        """
        return op.einsum('fF, nFx -> nxf', mask, fk)

    def build(self, input_shape):
        super().build(input_shape)
        if self.num_replicas > 1:
            self.compute_observable = compute_dis_observable_many_replica
        else:
            # Currying the function so that the `Observable` does not need
            # to get modified
            def compute_dis_observable_one_replica_w_pc(pdf, padded_fk):
                return compute_dis_observable_one_replica(pdf, padded_fk, power_corrections = self.power_corrections)
            self.compute_observable = compute_dis_observable_one_replica_w_pc


def compute_dis_observable_many_replica(pdf, padded_fk):
    """
    Contract masked fk table with PDF.

    Parameters
    ----------
        pdf: list[tensor]
            list of pdf of shape (batch=1, replicas, xgrid, flavours)
        padded_fk: tensor
            masked fk table of shape (ndata, xgrid, flavours)

    Returns
    -------
        tensor
            observable of shape (batch=1, replicas, ndata)
    """
    return op.einsum('brxf, nxf -> brn', pdf[0], padded_fk)


def compute_dis_observable_one_replica(pdf, padded_fk, power_corrections = None):
    """
    Same operations as above but a specialized implementation that is more efficient for 1 replica,
    masking the PDF rather than the fk table.
    """
    if power_corrections is None:

      return op.tensor_product(pdf, padded_fk, axes=[(2, 3), (1, 2)])
    else:

      return op.tensor_product(pdf, padded_fk, axes=[(2, 3), (1, 2)]) + power_corrections
