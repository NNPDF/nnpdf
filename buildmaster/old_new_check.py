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
    For instance commit 5372da065
"""

from argparse import ArgumentParser
from dataclasses import dataclass
from functools import cached_property
import sys

from lhapdf import setVerbosity
import numpy as np
import pandas as pd
import yaml

from validphys.api import API
from validphys.convolution import central_predictions
from validphys.datafiles import path_commondata as old_cd_root

dataset_names_path = old_cd_root.with_name("new_commondata") / "dataset_names.yml"

setVerbosity(0)
pdf = "NNPDF40_nnlo_as_01180"
pdf_load = API.pdf(pdf=pdf)


class CheckFailed(Exception):
    pass


@dataclass
class DToCompare:
    old_name: str
    new_name: str
    variant: str = None
    theory_id: int = 717

    def __str__(self):
        return f"[old: {old_name} vs new: {new_name}]"

    @property
    def generic(self):
        return {"use_cuts": "internal", "theoryid": self.theory_id, "pdf": pdf}

    @property
    def dataset_input_old(self):
        return {"dataset": self.old_name}

    @property
    def dataset_input_new(self):
        if self.variant is None:
            return {"dataset": self.new_name}
        else:
            return {"dataset": self.new_name, "variant": self.variant}

    @cached_property
    def ds_old(self):
        """Load the old-commondata dataset"""
        return API.dataset(dataset_input=self.dataset_input_old, **self.generic)

    @cached_property
    def cd_old(self):
        """Load the commondata object out of the dataset"""
        return self.ds_old.load_commondata()

    @cached_property
    def ds_new(self):
        """Load the new-commondata dataset"""
        return API.dataset(dataset_input=self.dataset_input_new, **self.generic)

    @cached_property
    def cd_new(self):
        """Load the commondata object out of the dataset"""
        return self.ds_new.load_commondata()

    def api_call(self, api_function, extra_config=None):
        """Apply a validphys API call for a given function for both the old and new datasets"""
        if extra_config is None:
            extra_config = {}
        old_val = getattr(API, api_function)(
            dataset_input=self.dataset_input_old, **self.generic, **extra_config
        )
        new_val = getattr(API, api_function)(
            dataset_input=self.dataset_input_new, **self.generic, **extra_config
        )
        return old_val, new_val


def check_central_values(dcontainer, strict=True):
    """Check the central values

    By default, exit if they are not _the same_
    """
    cd_old = dcontainer.cd_old
    cd_new = dcontainer.cd_new
    if np.allclose(cd_old.central_values, cd_new.central_values):
        return True

    print(f"# Problem in the data comparison of {dcontainer}")
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


def check_chi2(dcontainer, strict=True, rtol=1e-5):
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
        raise CheckFailed

    print(f"# Differences in the computation of chi2 {chi2_old:.5} vs {chi2_new:.5}")
    # Check the predictions first
    new_pred = central_predictions(dcontainer.ds_new, pdf_load)
    old_pred = central_predictions(dcontainer.ds_old, pdf_load)
    if not np.allclose(new_pred, old_pred):
        print("... but the predictions were already different")

    old_covmat, new_covmat = dcontainer.api_call("covmat_from_systematics")
    if not np.allclose(old_covmat, new_covmat):
        print("    The covmats are different", end="")
        if not np.allclose(np.diag(old_covmat), np.diag(new_covmat)):
            print(" ...even the diagonal!")
        else:
            print(" ...but the diagonal is the same!")


def run_comparison(old_name, new_name, variant=None, theory_id=717):
    """
    Given an old and new datasets, perform a full comparison
    """
    # Create the DToCompare container class
    # which knows how to call validphys for the various informations it needs
    # and eases the printing
    dcontainer = DToCompare(old_name, new_name, variant=variant, theory_id=theory_id)

    # Check central data
    check_central_values(dcontainer)

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
        print(f"The t0 chi2 is different: {chi2_old_t0:.5} vs {chi2_new_t0:.5}")
        raise CheckFailed

    print(f"> Comparison ok! {dcontainer}")


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

        try:
            run_comparison(old_name, new_name, variant, theory_id=args.tid)
        except CheckFailed as e:
            if args.stop:
                raise e
        except Exception as e:
            print(f"Failure for {old_name}: {new_name}")
            if args.stop:
                raise e
            if args.verbose:
                print(e)
        except BaseException as e:
            print(f"Failure for {old_name}: {new_name}")
            if args.stop:
                raise e
            if args.verbose:
                print(e)
