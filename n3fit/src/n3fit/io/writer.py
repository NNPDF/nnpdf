"""
    Module containing functions dedicated to the write down of the output of n3fit

    The goal is to generate the same folder/file structure as the old nnfit code
    so previously active scripts can still work.
"""

import os
import time
import json
import numpy as np
from reportengine.compat import yaml
import validphys
from validphys.pdfgrids import EXPORT_XGRID
import n3fit


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
        # Compute the arclengths
        arc_lengths = self.pdf_object.compute_arclength()
        # Compute the integrability numbers
        integrability_numbers = self.pdf_object.integrability_numbers()
        # Construct the chi2exp file
        allchi2_lines = self.stopping_object.chi2exps_str()
        # Construct the preproc file
        preproc_lines = ""  # TODO decide how to fill up this in a sensible way

        # Check the directory exist, if it doesn't, generate it
        os.makedirs(replica_path_set, exist_ok=True)

        # export PDF grid to file
        storefit(
            self.pdf_object,
            self.replica_number,
            replica_path_set,
            fitname,
            self.q2,
            self.stopping_object.e_best_chi2,
            vl_chi2,
            tr_chi2,
            true_chi2,
            arc_lengths,
            integrability_numbers,
            allchi2_lines,
            preproc_lines,
            self.stopping_object.positivity_status(),
            self.timings,
        )

        # TODO: compute the chi2s directly from the stopping object
        # export all metadata from the fit to a single yaml file
        output_file = f"{replica_path_set}/{fitname}.json"
        json_dict = jsonfit(
            self.stopping_object, self.pdf_object, tr_chi2, vl_chi2, true_chi2, self.timings
        )
        with open(output_file, "w") as fs:
            json.dump(json_dict, fs, indent=2)


def jsonfit(stopping_object, pdf_object, tr_chi2, vl_chi2, true_chi2, timing):
    """Generates a dictionary containing all relevant metadata for the fit

    Parameters
    ----------
        stopping_object: n3fit.stopping.Stopping
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
        timing: dict
            dictionary of the timing of the different events that happened
    """
    all_info = {}
    # Generate preprocessing information
    all_info["preprocessing"] = pdf_object.get_preprocessing_factors()
    # .fitinfo-like info
    all_info["stop_epoch"] = stopping_object.stop_epoch
    all_info["best_epoch"] = stopping_object.e_best_chi2
    all_info["erf_tr"] = tr_chi2
    all_info["erf_vl"] = vl_chi2
    all_info["chi2"] = true_chi2
    all_info["pos_state"] = stopping_object.positivity_status()
    all_info["arc_lengths"] = pdf_object.compute_arclength().tolist()
    all_info["integrability"] = pdf_object.integrability_numbers().tolist()
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
    nite,
    erf_vl,
    erf_tr,
    chi2,
    arc_lengths,
    integrability_numbers,
    all_chi2_lines,
    all_preproc_lines,
    pos_state,
    timings,
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
        `nite`
            number of epochs before stopping
        `erf_vl`
            chi2 for the validation
        `erf_tr`
            chi2 for the training
        `chi2`
            chi2 for the exp (unreplica'd data)
        `arc_lengths`
            list with the arclengths of all replicas
        `integrability_numbers`
            list with the integrability numbers
        `all_chi2_lines`
            a big log with the chi2 every 100 epochs
        `all_preproc_lines`
            (None)
        `pos_state`
            positivity passing flag
        `timings`
            dictionary to dump to file containing timing information of the fit
    """
    result = pdf_object(EXPORT_XGRID, flavours="n3fit")
    lha = evln2lha(result.T).T

    data = {
        "replica": replica,
        "q20": q20,
        "xgrid": EXPORT_XGRID.T.tolist()[0],
        "labels": ["TBAR", "BBAR", "CBAR", "SBAR", "UBAR", "DBAR", "GLUON", "D", "U", "S", "C", "B", "T", "PHT"],
        "pdfgrid": lha.tolist(),
    }

    with open(f"{replica_path}/{fitname}.exportgrid", "w") as fs:
        yaml.dump(data, fs)

    # create empty files to make postfit happy
    emptyfiles = ["chi2exps.log", f"{fitname}.params", f"{fitname}.sumrules"]
    for fs in emptyfiles:
        open(f"{replica_path}/{fs}", "a").close()

    # Write chi2exp
    with open(f"{replica_path}/chi2exps.log", "w") as fs:
        for line in all_chi2_lines:
            fs.write(line)

    # Write preproc information
    with open(f"{replica_path}/{fitname}.preproc", "w") as fs:
        for line in all_preproc_lines:
            fs.write(line)

    # create info file
    arc_line = " ".join(str(i) for i in arc_lengths)
    integrability_line = " ".join(str(i) for i in integrability_numbers)
    with open(f"{replica_path}/{fitname}.fitinfo", "w") as fs:
        fs.write(f"{nite} {erf_vl} {erf_tr} {chi2} {pos_state}\n")
        fs.write(arc_line)
        fs.write(f"\n")
        fs.write(integrability_line)

    # create .time file
    with open(f"{replica_path}/{fitname}.time", "w") as fs:
        json.dump(timings, fs, indent=2)

    # create .version file, with the version of programs used
    with open(f"{replica_path}/version.info", "w") as fs:
        versions = version()
        json.dump(versions, fs, indent=2)
