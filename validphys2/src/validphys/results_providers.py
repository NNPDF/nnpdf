"""
results_providers.py

module which bridges between underlying functions concerned with loading
theory predictions and data and actions which can be accessed
by other actions/providers.

"""
from reportengine import collect

from validphys.commondataparser import load_commondata
from validphys.convolution import central_predictions

def loaded_commondata_with_cuts(commondata, cuts):
    """Load the commondata and apply cuts.

    Parameters
    ----------
    commondata: validphys.core.CommonDataSpec
        commondata to load and cut.
    cuts: validphys.core.cuts, None
        valid cuts, used to cut loaded commondata.

    Returns
    -------
    loaded_cut_commondata: validphys.coredata.CommonData

    """
    lcd = load_commondata(commondata)
    return lcd.with_cuts(cuts)

dataset_inputs_loaded_cd_with_cuts = collect(
    "loaded_commondata_with_cuts", ("data_input",))


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
