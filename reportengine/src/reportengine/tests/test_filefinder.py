#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 12:32:31 2017

@author: zah
"""
import pathlib
import shutil

import pytest

import reportengine
from reportengine.filefinder import FallbackFinder, ModuleFinder

from reportengine.tests.utils import tmp

def test_fileloader(tmp):
    p = tmp
    patata = (p/'patata')
    patata.touch()
    d = p/'dir'
    d.mkdir()
    (d/'file').touch()


    f = FallbackFinder([tmp])
    assert f.find('patata') == (p, 'patata')
    assert f.find('dir/file') == (p, 'dir/file')

    assert set(f.hint_files()) == {(p, 'patata'), (p, 'dir')}

    with pytest.raises(ValueError):
        FallbackFinder([patata])
    shutil.rmtree(p)
    with pytest.raises(FileNotFoundError):
        f.find('patata')

    mf = ModuleFinder(reportengine)
    assert mf.find('__init__.py')
