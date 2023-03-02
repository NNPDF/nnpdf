"""
    Module containing functions dedicated to the write down of the output of n3fit

    The goal is to generate the same folder/file structure as the old nnfit code
    so previously active scripts can still work.
"""
import os
import json
import numpy as np
from reportengine.compat import yaml
import validphys
import n3fit
from n3fit import vpinterface
from evolven3fit_new.cli import XGRID

class WriterWrapper:
    def __init__(self, replica_number, pdf_object, stopping_object, q2, timings):
        """
        Initializes the writer for one given replica. This is decoupled from the writing
        of the fit in order to fix some of the variables which would be, in principle,
        be shared by several different history objects.

        Parameters
        ----------
            `replica_number`
                index of the replica
            `pdf_object`
                function to evaluate with a grid in x to generate a pdf
            `stopping_object`
                a stopping.Stopping object
            `q2`
                q^2 of the fit
            `timings`
                dictionary of the timing of the different events that happened
        """
        self.replica_number = replica_number
        self.pdf_object = pdf_object
        self.stopping_object = stopping_object
        self.q2 = q2
        self.timings = timings

    def write_data(self, replica_path_set, fitname, tr_chi2, vl_chi2, true_chi2):
        """
        Wrapper around the `storefit` function.

        Parameters
        ----------
            `replica_path_set`
                full path for the replica, ex: `${PWD}/runcard_name/nnfit/replica_1`
            `fitname`
                name of the fit
            `tr_chi2`
                training chi2
            `vl_chi2`
                validation chi2
            `true_chi2`
                chi2 of the replica to the central experimental data
        """
        # Check the directory exist, if it doesn't, generate it
        os.makedirs(replica_path_set, exist_ok=True)

        stop_epoch = self.stopping_object.stop_epoch

        # Get the replica status for this object
        replica_status = self.stopping_object.get_next_replica()

        # export PDF grid to file
        storefit(
            self.pdf_object,
            self.replica_number,
            replica_path_set,
            fitname,
            self.q2,
        )

        # write the log file for the chi2
        chi2_log = self.stopping_object.chi2exps_json()
        with (replica_path_set / "chi2exps.log").open("w", encoding="utf-8") as fs:
            json.dump(chi2_log, fs, indent=2, cls = SuperEncoder)

        # export all metadata from the fit to a single yaml file
        output_file = f"{replica_path_set}/{fitname}.json"
        json_dict = jsonfit(
            replica_status, self.pdf_object, tr_chi2, vl_chi2, true_chi2, stop_epoch, self.timings
        )
        with open(output_file, "w", encoding="utf-8") as fs:
            json.dump(json_dict, fs, indent=2, cls = SuperEncoder)


class SuperEncoder(json.JSONEncoder):
    """ Custom json encoder to get around the fact that np.float32 =/= float """
    def default(self, o):
        if isinstance(o, np.float32):
            return float(o)
        return super().default(o)


def jsonfit(replica_status, pdf_object, tr_chi2, vl_chi2, true_chi2, stop_epoch, timing):
    """Generates a dictionary containing all relevant metadata for the fit

    Parameters
    ----------
        replica_status: n3fit.stopping.ReplicaBest
            a stopping.Validation object
        pdf_object: n3fit.vpinterface.N3PDF
            N3PDF object constructed from the pdf_model
            that receives as input a point in x and returns an array of 14 flavours
        tr_chi2: float
            chi2 for the training
        vl_chi2: float
            chi2 for the validation
        true_chi2: float
            chi2 for the exp (unreplica'd data)
        epoch_stop: int
            epoch at which the stopping stopped (not the one for the best fit!)
        timing: dict
            dictionary of the timing of the different events that happened
    """
    all_info = {}
    # Generate preprocessing information
    all_info["preprocessing"] = pdf_object.get_preprocessing_factors()
    # .fitinfo-like info
    all_info["stop_epoch"] = stop_epoch
    all_info["best_epoch"] = replica_status.best_epoch
    all_info["erf_tr"] = tr_chi2
    all_info["erf_vl"] = vl_chi2
    all_info["chi2"] = true_chi2
    all_info["pos_state"] = replica_status.positivity_status
    all_info["arc_lengths"] = vpinterface.compute_arclength(pdf_object).tolist()
    all_info["integrability"] = vpinterface.integrability_numbers(pdf_object).tolist()
    all_info["timing"] = timing
    # Versioning info
    all_info["version"] = version()
    return all_info


def version():
    """ Generates a dictionary with misc version info for this run """
    versions = {}
    try:
        # Wrap tf in try-except block as it could possible to run n3fit without tf
        import tensorflow as tf
        from tensorflow.python.framework import test_util

        versions["keras"] = tf.keras.__version__
        mkl = test_util.IsMklEnabled()
        versions["tensorflow"] = f"{tf.__version__}, mkl={mkl}"
    except ImportError:
        versions["tensorflow"] = "Not available"
        versions["keras"] = "Not available"
    except AttributeError:
        # Check for MKL was only recently introduced and is not part of the official API
        versions["tensorflow"] = f"{tf.__version__}, mkl=??"
    except:
        # We don't want _any_ uncaught exception to crash the whole program at this point
        pass
    versions["numpy"] = np.__version__
    versions["nnpdf"] = n3fit.__version__
    try:
        versions["validphys"] = validphys.__version__
    except AttributeError:
        versions["validphys"] = "unknown"
    return versions


def evln2lha(evln):
    # evln Basis
    # {"PHT","SNG","GLU","VAL","V03","V08","V15","V24","V35","T03","T08","T15","T24","T35"};
    # lha Basis:
    # {"TBAR","BBAR","CBAR","SBAR","UBAR","DBAR","GLUON","D","U","S","C","B","T","PHT"}
    lha = np.zeros(evln.shape)
    lha[13] = evln[0]

    lha[6] = evln[2]

    lha[8] = ( 10*evln[1]
	       + 30*evln[9] + 10*evln[10] + 5*evln[11] + 3*evln[12] + 2*evln[13]
	       + 10*evln[3] + 30*evln[4] + 10*evln[5] + 5*evln[6] + 3*evln[7] + 2*evln[8] ) / 120

    lha[4] = ( 10*evln[1]
		  + 30*evln[9] + 10*evln[10] + 5*evln[11] + 3*evln[12] + 2*evln[13]
		  - 10*evln[3] - 30*evln[4] - 10*evln[5] - 5*evln[6] - 3*evln[7] - 2*evln[8] ) / 120

    lha[7] = ( 10*evln[1]
	       - 30*evln[9] + 10*evln[10] + 5*evln[11] + 3*evln[12] + 2*evln[13]
	       + 10*evln[3] - 30*evln[4] + 10*evln[5] + 5*evln[6] + 3*evln[7] + 2*evln[8] ) / 120

    lha[5] = ( 10*evln[1]
		  - 30*evln[9] + 10*evln[10] + 5*evln[11] + 3*evln[12] + 2*evln[13]
		  - 10*evln[3] + 30*evln[4] - 10*evln[5] - 5*evln[6] - 3*evln[7] - 2*evln[8] ) / 120

    lha[9] = ( 10*evln[1]
	       - 20*evln[10] + 5*evln[11] + 3*evln[12] + 2*evln[13]
	       + 10*evln[3] - 20*evln[5] + 5*evln[6] + 3*evln[7] + 2*evln[8] ) / 120

    lha[3] = ( 10*evln[1]
		  - 20*evln[10] + 5*evln[11] + 3*evln[12] + 2*evln[13]
		  - 10*evln[3] + 20*evln[5] - 5*evln[6] - 3*evln[7] - 2*evln[8] ) / 120

    lha[10] = ( 10*evln[1]
	       - 15*evln[11] + 3*evln[12] + 2*evln[13]
	       + 10*evln[3] - 15*evln[6] + 3*evln[7] + 2*evln[8] ) / 120

    lha[2] = ( 10*evln[1]
		  - 15*evln[11] + 3*evln[12] + 2*evln[13]
		  - 10*evln[3] + 15*evln[6] - 3*evln[7] - 2*evln[8] ) / 120

    lha[11] = ( 5*evln[1]
	       - 6*evln[12] + evln[13]
	       + 5*evln[3] - 6*evln[7] + evln[8] ) / 60

    lha[1] = ( 5*evln[1]
		  - 6*evln[12] + evln[13]
		  - 5*evln[3] + 6*evln[7] - evln[8] ) / 60

    lha[12] = ( evln[1]
	       - evln[13]
	       + evln[3] - evln[8] ) / 12

    lha[0] = ( evln[1]
		  - evln[13]
		  - evln[3] + evln[8] ) / 12
    return lha


def storefit(
    pdf_object,
    replica,
    replica_path,
    fitname,
    q20,
):
    """
    One-trick function which generates all output in the NNPDF format
    so that all other scripts can still be used.

    Parameters
    ----------
        `pdf_object`
            N3PDF object constructed from the pdf_model
            that receives as input a point in x and returns an array of 14 flavours
        `replica`
            the replica index
        `replica_path`
            path for this replica
        `fitname`
            name of the fit
        `q20`
            q_0^2
    """
    # build exportgrid
    xgrid = XGRID.reshape(-1, 1)
        
    result = pdf_object(xgrid, flavours="n3fit").squeeze()
    lha = evln2lha(result.T).T

    data = {
        "replica": replica,
        "q20": q20,
        "xgrid": xgrid.tolist(),
        "labels": ["TBAR", "BBAR", "CBAR", "SBAR", "UBAR", "DBAR", "GLUON", "D", "U", "S", "C", "B", "T", "PHT"],
        "pdfgrid": lha.tolist(),
    }

    with open(f"{replica_path}/{fitname}.exportgrid", "w") as fs:
        yaml.dump(data, fs)
