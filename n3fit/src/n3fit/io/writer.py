"""
    Module containing functions dedicated to the write down of the output of n3fit

    The goal is to generate the same folder/file structure as the old nnfit code
    so previously active scripts can still work.
"""
import json
import logging

import numpy as np

import n3fit
from n3fit import vpinterface
from reportengine.compat import yaml
import validphys

log = logging.getLogger(__name__)

XGRID = np.array(
    [
        1.00000000000000e-09,
        1.29708482343957e-09,
        1.68242903474257e-09,
        2.18225315420583e-09,
        2.83056741739819e-09,
        3.67148597892941e-09,
        4.76222862935315e-09,
        6.17701427376180e-09,
        8.01211109898438e-09,
        1.03923870607245e-08,
        1.34798064073805e-08,
        1.74844503691778e-08,
        2.26788118881103e-08,
        2.94163370300835e-08,
        3.81554746595878e-08,
        4.94908707232129e-08,
        6.41938295708371e-08,
        8.32647951986859e-08,
        1.08001422993829e-07,
        1.40086873081130e-07,
        1.81704331793772e-07,
        2.35685551545377e-07,
        3.05703512595323e-07,
        3.96522309841747e-07,
        5.14321257236570e-07,
        6.67115245136676e-07,
        8.65299922973143e-07,
        1.12235875241487e-06,
        1.45577995547683e-06,
        1.88824560514613e-06,
        2.44917352454946e-06,
        3.17671650028717e-06,
        4.12035415232797e-06,
        5.34425265752090e-06,
        6.93161897806315e-06,
        8.99034258238145e-06,
        1.16603030112258e-05,
        1.51228312288769e-05,
        1.96129529349212e-05,
        2.54352207134502e-05,
        3.29841683435992e-05,
        4.27707053972016e-05,
        5.54561248105849e-05,
        7.18958313632514e-05,
        9.31954227979614e-05,
        1.20782367731330e-04,
        1.56497209466554e-04,
        2.02708936328495e-04,
        2.62459799331951e-04,
        3.39645244168985e-04,
        4.39234443000422e-04,
        5.67535660104533e-04,
        7.32507615725537e-04,
        9.44112105452451e-04,
        1.21469317686978e-03,
        1.55935306118224e-03,
        1.99627451141338e-03,
        2.54691493736552e-03,
        3.23597510213126e-03,
        4.09103436509565e-03,
        5.14175977083962e-03,
        6.41865096062317e-03,
        7.95137940306351e-03,
        9.76689999624100e-03,
        1.18876139251364e-02,
        1.43298947643919e-02,
        1.71032279460271e-02,
        2.02100733925079e-02,
        2.36463971369542e-02,
        2.74026915728357e-02,
        3.14652506132444e-02,
        3.58174829282429e-02,
        4.04411060163317e-02,
        4.53171343973807e-02,
        5.04266347950069e-02,
        5.57512610084339e-02,
        6.12736019390519e-02,
        6.69773829498255e-02,
        7.28475589986517e-02,
        7.88703322292727e-02,
        8.50331197801452e-02,
        9.13244910278679e-02,
        9.77340879783772e-02,
        1.04252538208639e-01,
        1.10871366547237e-01,
        1.17582909372878e-01,
        1.24380233801599e-01,
        1.31257062945031e-01,
        1.38207707707289e-01,
        1.45227005135651e-01,
        1.52310263065985e-01,
        1.59453210652156e-01,
        1.66651954293987e-01,
        1.73902938455578e-01,
        1.81202910873333e-01,
        1.88548891679097e-01,
        1.95938145999193e-01,
        2.03368159629765e-01,
        2.10836617429103e-01,
        2.18341384106561e-01,
        2.25880487124065e-01,
        2.33452101459503e-01,
        2.41054536011681e-01,
        2.48686221452762e-01,
        2.56345699358723e-01,
        2.64031612468684e-01,
        2.71742695942783e-01,
        2.79477769504149e-01,
        2.87235730364833e-01,
        2.95015546847664e-01,
        3.02816252626866e-01,
        3.10636941519503e-01,
        3.18476762768082e-01,
        3.26334916761672e-01,
        3.34210651149156e-01,
        3.42103257303627e-01,
        3.50012067101685e-01,
        3.57936449985571e-01,
        3.65875810279643e-01,
        3.73829584735962e-01,
        3.81797240286494e-01,
        3.89778271981947e-01,
        3.97772201099286e-01,
        4.05778573402340e-01,
        4.13796957540671e-01,
        4.21826943574548e-01,
        4.29868141614175e-01,
        4.37920180563205e-01,
        4.45982706956990e-01,
        4.54055383887562e-01,
        4.62137890007651e-01,
        4.70229918607142e-01,
        4.78331176755675e-01,
        4.86441384506059e-01,
        4.94560274153348e-01,
        5.02687589545177e-01,
        5.10823085439086e-01,
        5.18966526903235e-01,
        5.27117688756998e-01,
        5.35276355048428e-01,
        5.43442318565661e-01,
        5.51615380379768e-01,
        5.59795349416641e-01,
        5.67982042055800e-01,
        5.76175281754088e-01,
        5.84374898692498e-01,
        5.92580729444440e-01,
        6.00792616663950e-01,
        6.09010408792398e-01,
        6.17233959782450e-01,
        6.25463128838069e-01,
        6.33697780169485e-01,
        6.41937782762089e-01,
        6.50183010158361e-01,
        6.58433340251944e-01,
        6.66688655093089e-01,
        6.74948840704708e-01,
        6.83213786908386e-01,
        6.91483387159697e-01,
        6.99757538392251e-01,
        7.08036140869916e-01,
        7.16319098046733e-01,
        7.24606316434025e-01,
        7.32897705474271e-01,
        7.41193177421404e-01,
        7.49492647227008e-01,
        7.57796032432224e-01,
        7.66103253064927e-01,
        7.74414231541921e-01,
        7.82728892575836e-01,
        7.91047163086478e-01,
        7.99368972116378e-01,
        8.07694250750291e-01,
        8.16022932038457e-01,
        8.24354950923382e-01,
        8.32690244169987e-01,
        8.41028750298844e-01,
        8.49370409522600e-01,
        8.57715163684985e-01,
        8.66062956202683e-01,
        8.74413732009721e-01,
        8.82767437504206e-01,
        8.91124020497459e-01,
        8.99483430165226e-01,
        9.07845617001021e-01,
        9.16210532771399e-01,
        9.24578130473112e-01,
        9.32948364292029e-01,
        9.41321189563734e-01,
        9.49696562735755e-01,
        9.58074441331298e-01,
        9.66454783914439e-01,
        9.74837550056705e-01,
        9.83222700304978e-01,
        9.91610196150662e-01,
        1.00000000000000e00,
    ]
)


class WriterWrapper:
    def __init__(self, replica_numbers, pdf_objects, stopping_object, all_chi2s, theory, timings):
        """
        Initializes the writer for all replicas.

        This is decoupled from the writing of the fit in order to fix some of the variables
        which would be, in principle, be shared by several different history objects.

        Parameters
        ----------
            `replica_numbers`
                indices of the replicas
            `pdf_objects`
                function to evaluate with a grid in x to generate a pdf
            `stopping_object`
                a stopping.Stopping object
            `all_chi2s`
                list of all the chi2s, in the order: tr_chi2, vl_chi2, true_chi2
            `theory`
                theory information of the fit
            `timings`
                dictionary of the timing of the different events that happened
        """
        self.replica_numbers = replica_numbers
        self.pdf_objects = pdf_objects
        self.stopping_object = stopping_object
        self.theory = theory
        self.timings = timings
        self.tr_chi2, self.vl_chi2, self.true_chi2 = all_chi2s

    def write_data(self, save_path, fitname, weights_name):
        """
        Save all the data of a fit, for all replicas.

        Parameters
        ----------
            `save_path`
                path for the fit results, ex: `${PWD}/runcard_name/nnfit`
            `fitname`
                name of the fit, ex: `Basic_runcard`
            `weights_name`
                name of the file to save weights to, if not empty
        """
        save_path.mkdir(exist_ok=True, parents=True)

        self.preprocessing = []
        self.arc_lengths = []
        self.integrability_numbers = []
        for pdf_object in self.pdf_objects:
            self.preprocessing.append(pdf_object.get_preprocessing_factors())
            self.arc_lengths.append(vpinterface.compute_arclength(pdf_object).tolist())
            self.integrability_numbers.append(
                vpinterface.integrability_numbers(pdf_object).tolist()
            )

        for i, rn in enumerate(self.replica_numbers):
            replica_path = save_path / f"replica_{rn}"
            replica_path.mkdir(exist_ok=True, parents=True)

            self._write_chi2s(replica_path / "chi2exps.log")
            self._write_metadata_json(i, replica_path / f"{fitname}.json")
            self._export_pdf_grid(i, replica_path / f"{fitname}.exportgrid")
            if weights_name:
                self._write_weights(i, replica_path / f"{weights_name}")

    def _write_chi2s(self, out_path):
        # Note: same for all replicas, unless run separately
        chi2_log = self.stopping_object.chi2exps_json()
        with open(out_path, "w", encoding="utf-8") as fs:
            json.dump(chi2_log, fs, indent=2, cls=SuperEncoder)

    def _write_metadata_json(self, i, out_path):
        json_dict = jsonfit(
            best_epoch=self.stopping_object.e_best_chi2[i],
            positivity_status=self.stopping_object.positivity_statuses[i],
            preprocessing=self.preprocessing[i],
            arc_lengths=self.arc_lengths[i],
            integrability_numbers=self.integrability_numbers[i],
            tr_chi2=self.tr_chi2[i],
            vl_chi2=self.vl_chi2[i],
            true_chi2=self.true_chi2[i],
            # Note: the 2 arguments below are the same for all replicas, unless run separately
            timing=self.timings,
            stop_epoch=self.stopping_object.stop_epoch,
        )

        with open(out_path, "w", encoding="utf-8") as fs:
            json.dump(json_dict, fs, indent=2, cls=SuperEncoder)

        log.info(
            "Best fit for replica #%d, chi2=%.3f (tr=%.3f, vl=%.3f), at epoch %d.",
            self.replica_numbers[i],
            self.true_chi2[i],
            self.tr_chi2[i],
            self.vl_chi2[i],
            self.stopping_object.e_best_chi2[i],
        )

    def _export_pdf_grid(self, i, out_path):
        storefit(self.pdf_objects[i], self.replica_numbers[i], out_path, self.theory)

    def _write_weights(self, i, out_path):
        log.info(" > Saving the weights for future in %s", out_path)
        # Extract model out of N3PDF
        model = self.pdf_objects[i]._models[0]
        model.save_weights(out_path, save_format="h5")


class SuperEncoder(json.JSONEncoder):
    """Custom json encoder to get around the fact that np.float32 =/= float"""

    def default(self, o):
        if isinstance(o, np.float32):
            return float(o)
        return super().default(o)


def jsonfit(
    best_epoch,
    positivity_status,
    preprocessing,
    arc_lengths,
    integrability_numbers,
    tr_chi2,
    vl_chi2,
    true_chi2,
    stop_epoch,
    timing,
):
    """Generates a dictionary containing all relevant metadata for the fit

    Parameters
    ----------
        best_epoch: int
            epoch at which the best fit was found
        positivity_status: str
            string describing the positivity status of the fit
        preprocessing: dict
            dictionary of the preprocessing factors
        arc_lengths: list
            list of the arc lengths of the different PDFs
        integrability_numbers: list
            list of the integrability numbers of the different PDFs
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
    all_info["preprocessing"] = preprocessing
    # .fitinfo-like info
    all_info["stop_epoch"] = stop_epoch
    all_info["best_epoch"] = best_epoch
    all_info["erf_tr"] = tr_chi2
    all_info["erf_vl"] = vl_chi2
    all_info["chi2"] = true_chi2
    all_info["pos_state"] = positivity_status
    all_info["arc_lengths"] = arc_lengths
    all_info["integrability"] = integrability_numbers
    all_info["timing"] = timing
    # Versioning info
    all_info["version"] = version()
    return all_info


def version():
    """Generates a dictionary with misc version info for this run"""
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


def evln2lha(evln, nf=6):
    # evln Basis
    # {"PHT","SNG","GLU","VAL","V03","V08","V15","V24","V35","T03","T08","T15","T24","T35"};
    # lha Basis:
    # {"TBAR","BBAR","CBAR","SBAR","UBAR","DBAR","GLUON","D","U","S","C","B","T","PHT"}
    lha = np.zeros(evln.shape)
    lha[13] = evln[0]

    lha[6] = evln[2]

    lha[8] = (
        10 * evln[1]
        + 30 * evln[9]
        + 10 * evln[10]
        + 5 * evln[11]
        + 3 * evln[12]
        + 2 * evln[13]
        + 10 * evln[3]
        + 30 * evln[4]
        + 10 * evln[5]
        + 5 * evln[6]
        + 3 * evln[7]
        + 2 * evln[8]
    ) / 120

    lha[4] = (
        10 * evln[1]
        + 30 * evln[9]
        + 10 * evln[10]
        + 5 * evln[11]
        + 3 * evln[12]
        + 2 * evln[13]
        - 10 * evln[3]
        - 30 * evln[4]
        - 10 * evln[5]
        - 5 * evln[6]
        - 3 * evln[7]
        - 2 * evln[8]
    ) / 120

    lha[7] = (
        10 * evln[1]
        - 30 * evln[9]
        + 10 * evln[10]
        + 5 * evln[11]
        + 3 * evln[12]
        + 2 * evln[13]
        + 10 * evln[3]
        - 30 * evln[4]
        + 10 * evln[5]
        + 5 * evln[6]
        + 3 * evln[7]
        + 2 * evln[8]
    ) / 120

    lha[5] = (
        10 * evln[1]
        - 30 * evln[9]
        + 10 * evln[10]
        + 5 * evln[11]
        + 3 * evln[12]
        + 2 * evln[13]
        - 10 * evln[3]
        + 30 * evln[4]
        - 10 * evln[5]
        - 5 * evln[6]
        - 3 * evln[7]
        - 2 * evln[8]
    ) / 120

    lha[9] = (
        10 * evln[1]
        - 20 * evln[10]
        + 5 * evln[11]
        + 3 * evln[12]
        + 2 * evln[13]
        + 10 * evln[3]
        - 20 * evln[5]
        + 5 * evln[6]
        + 3 * evln[7]
        + 2 * evln[8]
    ) / 120

    lha[3] = (
        10 * evln[1]
        - 20 * evln[10]
        + 5 * evln[11]
        + 3 * evln[12]
        + 2 * evln[13]
        - 10 * evln[3]
        + 20 * evln[5]
        - 5 * evln[6]
        - 3 * evln[7]
        - 2 * evln[8]
    ) / 120

    # if a heavy quark is not active at Q0 (the scale at which the output of the fit is stored),
    # keep the PDF values at 0.0 to prevent small negative values due to numerical instabilities

    # charm
    if nf > 3:
        lha[10] = (
            10 * evln[1]
            - 15 * evln[11]
            + 3 * evln[12]
            + 2 * evln[13]
            + 10 * evln[3]
            - 15 * evln[6]
            + 3 * evln[7]
            + 2 * evln[8]
        ) / 120

        lha[2] = (
            10 * evln[1]
            - 15 * evln[11]
            + 3 * evln[12]
            + 2 * evln[13]
            - 10 * evln[3]
            + 15 * evln[6]
            - 3 * evln[7]
            - 2 * evln[8]
        ) / 120

    # bottom
    if nf > 4:
        lha[11] = (5 * evln[1] - 6 * evln[12] + evln[13] + 5 * evln[3] - 6 * evln[7] + evln[8]) / 60
        lha[1] = (5 * evln[1] - 6 * evln[12] + evln[13] - 5 * evln[3] + 6 * evln[7] - evln[8]) / 60

    # top
    if nf > 5:
        lha[12] = (evln[1] - evln[13] + evln[3] - evln[8]) / 12
        lha[0] = (evln[1] - evln[13] - evln[3] + evln[8]) / 12

    return lha


def storefit(pdf_object, replica, out_path, theory):
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
        `out_path`
            the path where to store the output
        `theory`
            theory information of the fit
    """
    # build exportgrid
    xgrid = XGRID.reshape(-1, 1)

    q0 = theory.get_description().get("Q0")
    q20 = q0**2

    # Set the active flavours, by default u, d, s
    active_flavours = 3

    # Now check whether c, b, t are active at Q0
    for quark in ["c", "b", "t"]:
        mass = theory.get_description().get(f"m{quark}")
        threshold = theory.get_description().get(f"k{quark}Thr")

        if q0 < mass * threshold:
            break
        active_flavours += 1

    # Double check IC
    if theory.get_description().get("IC") == 1:
        active_flavours = np.maximum(4, active_flavours)

    result = pdf_object(xgrid, flavours="n3fit").squeeze()
    lha = evln2lha(result.T, nf=active_flavours).T

    data = {
        "replica": replica,
        "q20": q20,
        "xgrid": xgrid.T.tolist()[0],
        "labels": [
            "TBAR",
            "BBAR",
            "CBAR",
            "SBAR",
            "UBAR",
            "DBAR",
            "GLUON",
            "D",
            "U",
            "S",
            "C",
            "B",
            "T",
            "PHT",
        ],
        "pdfgrid": lha.tolist(),
    }

    with open(out_path, "w") as fs:
        yaml.dump(data, fs)
