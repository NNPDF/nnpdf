from typing import Callable, List, Optional

import numpy as np
import numpy.typing as npt
from scipy.interpolate import PchipInterpolator


def generate_scaler(
    input_list: List[npt.NDArray], interpolation_points: Optional[int] = None
) -> Callable:
    """
    Generate the scaler function that applies feature scaling to the input data.

    Parameters
    ----------
    input_list : list of numpy.ndarray
        The list of input data arrays.
    interpolation_points : int, optional

    Returns
    -------
    _scaler : Callable
        The scaler function that applies feature scaling to the input data.
    """
    input_arr = np.concatenate(input_list, axis=1)
    input_arr = np.sort(input_arr)
    input_arr_size = input_arr.size

    # Define an evenly spaced grid in the domain [0,1]
    # force_set_smallest is used to make sure the smallest point included in the scaling is
    # 1e-9, to prevent trouble when saving it to the LHAPDF grid
    force_set_smallest = input_arr.min() > 1e-9
    if force_set_smallest:
        new_xgrid = np.linspace(
            start=1 / input_arr_size, stop=1.0, endpoint=False, num=input_arr_size
        )
    else:
        new_xgrid = np.linspace(start=0, stop=1.0, endpoint=False, num=input_arr_size)

    # When mapping the FK xgrids onto our new grid, we need to consider degeneracies among
    # the x-values in the FK grids
    unique, counts = np.unique(input_arr, return_counts=True)
    map_to_complete = []
    for cumsum_ in np.cumsum(counts):
        # Make sure to include the smallest new_xgrid value, such that we have a point at
        # x<=1e-9
        map_to_complete.append(new_xgrid[cumsum_ - counts[0]])
    map_to_complete = np.array(map_to_complete)
    map_from_complete = unique

    #  If needed, set feature_scaling(x=1e-9)=0
    if force_set_smallest:
        map_from_complete = np.insert(map_from_complete, 0, 1e-9)
        map_to_complete = np.insert(map_to_complete, 0, 0.0)

    # Select the indices of the points that will be used by the interpolator
    onein = map_from_complete.size / (int(interpolation_points) - 1)
    selected_points = [round(i * onein - 1) for i in range(1, int(interpolation_points))]
    if selected_points[0] != 0:
        selected_points = [0] + selected_points
    map_from = map_from_complete[selected_points]
    map_from = np.log(map_from)
    map_to = map_to_complete[selected_points]

    try:
        scaler = PchipInterpolator(map_from, map_to)
    except ValueError as e:
        raise ValueError(
            "interpolation_points is larger than the number of unique input x-values"
        ) from e

    def _scaler(x):
        x_scaled = scaler(np.log(x))
        x_scaled = 2 * x_scaled - 1
        return np.concatenate([x_scaled, x], axis=-1)

    return _scaler
