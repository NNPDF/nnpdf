"""Tests of the QED couplings of the theory cards"""

import dataclasses

import pytest

from nnpdf_data import THEORY_CARDS_PATH
from nnpdf_data.theory import TheoryCard, TheoryCardError, qed_running_alpha
from nnpdf_data.theorydbutils import fetch_theory

# NNPDF4.1 NNLO + QED: G_mu scheme, alphaqed = alpha_Gmu = 1/132.10
THEORY_GMU = 41_041_001


def _card(**changes):
    theory = dict(fetch_theory(THEORY_CARDS_PATH, THEORY_GMU))  # do not modify the cached card
    theory.update(changes)
    fields = {f.name for f in dataclasses.fields(TheoryCard)}
    return TheoryCard(**{k: v for k, v in theory.items() if k in fields})


def test_alpha_qed_running_mz():
    """With alpha_qed_running_mz the QED running uses it, otherwise alphaqed (as in NNPDF4.0)"""
    card = _card(alpha_qed_running_mz=1 / 128.94)
    assert qed_running_alpha(card.asdict()) == 1 / 128.94
    card = _card()
    assert qed_running_alpha(card.asdict()) == card.alphaqed


def test_alphaqed_is_gmu():
    """When alpha_qed_running_mz is given, alphaqed must follow from GF, MW and MZ"""
    _card(alpha_qed_running_mz=1 / 128.94)
    # the running coupling in place of alpha_Gmu is rejected
    with pytest.raises(TheoryCardError):
        _card(alpha_qed_running_mz=1 / 128.94, alphaqed=1 / 128.94)
    # as is an alpha_Gmu inconsistent with the masses (MW from PDG2022 instead of PDG2024)
    with pytest.raises(TheoryCardError):
        _card(alpha_qed_running_mz=1 / 128.94, MW=80.377)


def test_alpha_qed_running_mz_needs_qref_mz():
    """alpha_qed_running_mz is given at M_Z, so a card with Qref != MZ is rejected"""
    with pytest.raises(TheoryCardError):
        _card(alpha_qed_running_mz=1 / 128.94, Qref=100.0)
