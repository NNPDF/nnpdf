#!/usr/bin/env python3
"""
    This scripts checks that the datasets in the file ``dataset_names.yml``
    are truly the same the same for all intents and purposes
    For that we check that the:
        1. Central value are the same
        2. The covariance matrix are the same
        3. The computed chi2 is the same
        4. The t0 chi2 is the same (so multiplicative and additive are equivalent in the new and old implementation)

    In order to perform this check we require a version of validphys that has both the old and new leader.
    For instance commit 9ddf580f98c99f3c953e7be328df4a4c95c43c14

    All checks are done assuming the same theory for the old and new dataset.
"""

from argparse import ArgumentParser
from dataclasses import dataclass
from functools import cached_property

from lhapdf import setVerbosity
import numpy as np
import pandas as pd
import yaml

from validphys.api import API
from validphys.convolution import central_predictions
from nnpdf_data import path_commondata as old_cd_root
from validphys.loader import Loader

dataset_names_path = old_cd_root.with_name("new_commondata") / "dataset_names.yml"

setVerbosity(0)
pdf = "NNPDF40_nnlo_as_01180"
pdf_load = API.pdf(pdf=pdf)
DEFAULT_STRICT = False

# This will only work if the code was installed in edit mode, but this is a development script so, fair game
runcard40 = (
    old_cd_root.parent.parent.parent.parent.parent
    / "n3fit"
    / "runcards"
    / "reproduce_nnpdf40"
    / "NNPDF40_nnlo_as_01180_1000.yml"
)
if runcard40.exists():
    runcard40_data_raw = yaml.safe_load(runcard40.read_text())["dataset_inputs"]
    runcard40_data = [i["dataset"] for i in runcard40_data_raw]
else:
    runcard40_data = None


class CheckFailed(Exception):
    pass


@dataclass
class DToCompare:
    old_name: str
    new_name: str
    variant: str = None
    theory_id: int = 717
    old_theory_id: int = 717

    def __str__(self):
        return f"[old: {old_name} vs new: {new_name}]"

    @property
    def is_positivity(self):
        return self.new_name.startswith("NNPDF_POS")

    @property
    def is_integrability(self):
        return self.new_name.startswith("NNPDF_INTEG")

    @property
    def is_lagrange(self):
        return self.is_integrability or self.is_positivity

    @property
    def generic(self):
        return {"use_cuts": "internal", "pdf": pdf}

    @property
    def dataset_input_old(self):
        di = {"dataset": self.old_name}
        if self.is_lagrange:
            di["maxlambda"] = 1.0

        if self.is_positivity:
            return {"posdataset": di, "theoryid": self.old_theory_id}
        elif self.is_integrability:
            return {"integdataset": di, "theoryid": self.old_theory_id}
        return {"dataset_input": di, "theoryid": self.old_theory_id}

    @property
    def dataset_input_new(self):
        di = {"dataset": self.new_name}
        if self.variant is not None:
            di["variant"] = self.variant
        if self.is_lagrange:
            di["maxlambda"] = 1.0

        if self.is_positivity:
            return {"posdataset": di, "theoryid": self.theory_id}
        elif self.is_integrability:
            return {"integdataset": di, "theoryid": self.theory_id}
        return {"dataset_input": di, "theoryid": self.theory_id}

    def api_load_dataset(self, dinput):
        """Load a dataset (positivity or not) with VP"""
        if self.is_lagrange:
            if self.new_name.startswith("NNPDF_POS"):
                return API.posdataset(**self.generic, **dinput)
            else:
                return API.integdataset(**self.generic, **dinput)
        else:
            return API.dataset(**self.generic, **dinput)

    @cached_property
    def ds_old(self):
        """Load the old-commondata dataset"""
        return self.api_load_dataset(self.dataset_input_old)

    @cached_property
    def cd_old(self):
        """Load the commondata object out of the dataset"""
        return self.ds_old.load_commondata()

    @cached_property
    def ds_new(self):
        """Load the new-commondata dataset"""
        return self.api_load_dataset(self.dataset_input_new)

    @cached_property
    def cd_new(self):
        """Load the commondata object out of the dataset"""
        return self.ds_new.load_commondata()

    def api_call(self, api_function, extra_config=None):
        """Apply a validphys API call for a given function for both the old and new datasets"""
        if extra_config is None:
            extra_config = {}
        old_val = getattr(API, api_function)(
            **self.generic, **extra_config, **self.dataset_input_old
        )
        new_val = getattr(API, api_function)(
            **self.generic, **extra_config, **self.dataset_input_new
        )
        return old_val, new_val


def check_central_values(dcontainer, strict=DEFAULT_STRICT):
    """Check the central values

    By default, exit if they are not _the same_
    """
    cd_old = dcontainer.cd_old
    cd_new = dcontainer.cd_new
    if np.allclose(cd_old.central_values, cd_new.central_values):
        return True

    print(f"# Problem in the data values comparison of {dcontainer}")
    if strict:
        raise CheckFailed

    _od = cd_old.central_values
    _nd = cd_new.central_values
    rat = np.abs(_od / _nd)
    if not np.allclose(rat, 1.0, rtol=1e-3):
        if not np.allclose(rat, 1.0, rtol=1e-2):
            print("Relative differences are above 1e-2! Panic!")
            df = pd.concat([_od, _nd, rat], axis=1)
            breakpoint()
            return False
        else:
            print("Relative differences between 1e-3 and 1e-2... acceptable...")
    else:
        print("Relative differences under 1e-3, continuing comparison...")


def check_theory(dcontainer):
    """Returns the old and new predictions"""
    new_pred = central_predictions(dcontainer.ds_new, pdf_load)
    old_pred = central_predictions(dcontainer.ds_old, pdf_load)
    return old_pred, new_pred


def check_chi2(dcontainer, strict=DEFAULT_STRICT, rtol=1e-5):
    """Checks whether the chi2 is the same
    A failure in the comparison of the chi2 can come from either:
        1. Data
        2. Covmat
        3. Theory
    Try to give as much information as possible before failure
    """
    chi2_old, chi2_new = dcontainer.api_call("central_chi2")
    if np.allclose(chi2_old, chi2_new, rtol=rtol):
        return True

    if strict:
        raise CheckFailed(f"Different chi2: {chi2_old:.4} vs {chi2_new:.4}")

    print(f"# Differences in the computation of chi2 {chi2_old:.5} vs {chi2_new:.5}")
    # Check the predictions first
    old_pred, new_pred = check_theory(dcontainer)
    if not np.allclose(new_pred, old_pred):
        print("... but the predictions were already different")

    old_covmat, new_covmat = dcontainer.api_call("covmat_from_systematics")
    if not np.allclose(old_covmat, new_covmat):
        print("    The covmats are different", end="")
        if not np.allclose(np.diag(old_covmat), np.diag(new_covmat)):
            print(" ...even the diagonal!")
        else:
            print(" ...but the diagonal is the same!")


def run_comparison(dcontainer):
    """
    Given an old and new datasets in a container, perform a full comparison
    """
    # Check central data
    check_central_values(dcontainer)

    if dcontainer.is_lagrange:
        pred_old, pred_new = check_theory(dcontainer)
        if np.allclose(pred_old, pred_new):
            print(f" > Comparison ok for positivity dataset {dcontainer}")
            return

    # chi2!
    # Computing the chi2 is checking:
    #   1. That the theory is loaded in the same way
    #   2. That the (experimental) covmat is created equal
    # Note that any problem in the data check above might break this
    check_chi2(dcontainer)

    # Check the chi2... with t0!
    # This checks that the ADD and MULT systematics are loaded in the same way

    chi2_old_t0, chi2_new_t0 = dcontainer.api_call(
        "central_chi2", extra_config={"use_t0": True, "t0pdfset": pdf}
    )
    if not np.allclose(chi2_old_t0, chi2_new_t0, rtol=1e-5):
        raise CheckFailed(f"The t0 chi2 is different: {chi2_old_t0:.5} vs {chi2_new_t0:.5}")

    print(f"> Comparison ok! {dcontainer}")


def check_40(old_name):
    if runcard40_data is not None and old_name not in runcard40_data:
        print(f"\033[92m  Dataset {old_name} was not part of 4.0\033[0m")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-s", "--stop", help="Stop on failure by raising the exception", action="store_true"
    )
    parser.add_argument(
        "-l", "--only-legacy", help="Check only those with variant: legacy", action="store_true"
    )
    parser.add_argument(
        "-v", "--verbose", help="Print the whole information on the failure", action="store_true"
    )
    parser.add_argument(
        "-f",
        "--filter",
        help="Simple filter to select a subset of datasets, applied on the new data",
        type=str,
        nargs='+',
    )
    parser.add_argument("-t", "--tid", help="Theory id, default 717", default=717)
    parser.add_argument("--old_tid", help="Old Theory id, default 717", default=717)
    args = parser.parse_args()

    all_ds_names = yaml.safe_load(dataset_names_path.read_text())

    for old_name, new_ds in all_ds_names.items():
        if isinstance(new_ds, str):
            new_name = new_ds
            variant = None
        else:
            new_name = new_ds["dataset"]
            variant = new_ds.get("variant")

        if args.only_legacy and variant != "legacy":
            continue

        if args.filter is not None:
            if not any(filter_word in new_name for filter_word in args.filter):
                continue

        if args.verbose:
            print(f"###########\n")

        try:
            # Create the DToCompare container class
            # which knows how to call validphys for the various informations it needs
            # and eases the printing
            dcontainer = DToCompare(
                old_name, new_name, variant=variant, theory_id=args.tid, old_theory_id=args.tid
            )
            run_comparison(dcontainer)
        except CheckFailed as e:
            print(f"> Failure for \033[91m\033[1m{old_name}: {new_name}\033[0m\033[0m (check)")
            # Regardless of the failure mode, tell the user whether this is a 4.0 dataset
            # But only in case of failure, otherwise why should we care
            check_40(old_name)
            if args.stop:
                raise e
            if args.verbose:
                print(e)
        except Exception as e:
            print(f"> Failure for \033[91m\033[1m{old_name}: {new_name}\033[0m\033[0m")
            check_40(old_name)
            if args.stop:
                raise e
            if args.verbose:
                print(e)
        except BaseException as e:
            print(f"> Failure for \033[91m\033[1m{old_name}: {new_name}\033[0m\033[0m")
            # Before raising the exception, check whether this is just a question of not having the right theory names
            theory_info = Loader().check_theoryID(args.tid)
            fk_path = dcontainer.ds_new.commondata.metadata.theory.fktables_to_paths(
                theory_info.path
            )
            if not fk_path[0][0].exists():
                print("Seems like the theory dictionary is not pointing to actual theories")
            if args.stop:
                raise e
            if args.verbose:
                print(e)
