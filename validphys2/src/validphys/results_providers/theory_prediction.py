"""
theory_prediction.py

Module containing actions which return theory predictions associated with
datasets.

"""
from reportengine import collect

from validphys.convolution import central_predictions

def dataset_t0_predictions(dataset, t0set):
    """Returns the t0 predictions for a ``dataset`` which are the predictions
    calculated using the central member of ``pdf``. Note that if ``pdf`` has
    errortype ``replicas``, and the dataset is a hadronic observable then the
    predictions of the central member are subtly different to the central
    value of the replica predictions.

    Parameters
    ----------
    dataset: validphys.core.DataSetSpec
        dataset for which to calculate t0 predictions
    t0set: validphys.core.PDF
        pdf used to calculate the predictions

    Returns
    -------
    t0_predictions: np.array
        1-D numpy array with predictions for each of the cut datapoints.

    """
    # Squeeze values since t0_pred is DataFrame with shape n_data * 1
    return central_predictions(dataset, t0set).to_numpy().squeeze()

dataset_inputs_t0_predictions = collect("dataset_t0_predictions", ("data",))
