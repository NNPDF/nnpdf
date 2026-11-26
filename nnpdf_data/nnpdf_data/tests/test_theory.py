"""Tests of the check that alphaqed is the running QED coupling, not alpha_Gmu"""

import dataclasses
import math

import pytest

from nnpdf_data import THEORY_CARDS_PATH
from nnpdf_data.theory import TheoryCard, TheoryCardError
from nnpdf_data.theorydbutils import fetch_theory

# NNPDF4.0 QED theory, alphaqed = 1/128.94 at Qref = 91.2 GeV
THEORY_QED = 702


def _card(**changes):
    theory = dict(fetch_theory(THEORY_CARDS_PATH, THEORY_QED))  # do not modify the cached card
    theory.update(changes)
    fields = {f.name for f in dataclasses.fields(TheoryCard)}
    return TheoryCard(**{k: v for k, v in theory.items() if k in fields})


def _alpha_gmu(card):
    return math.sqrt(2) * card.GF * card.MW**2 * (1 - card.MW**2 / card.MZ**2) / math.pi


def test_running_alphaqed_accepted():
    """The running coupling (NNPDF4.0 value and the MSbar PDG value) passes the check"""
    _card()
    _card(alphaqed=1 / 127.93)


def test_alpha_gmu_rejected():
    """With QED on, alphaqed equal to alpha_Gmu is rejected"""
    alpha_gmu = _alpha_gmu(_card())
    with pytest.raises(TheoryCardError):
        _card(alphaqed=alpha_gmu)
    with pytest.raises(TheoryCardError):
        _card(alphaqed=alpha_gmu * 1.01)


def test_alpha_gmu_allowed_without_qed():
    """Without QED, alphaqed is not used by the evolution, so it is not checked"""
    _card(QED=0, alphaqed=_alpha_gmu(_card()))
