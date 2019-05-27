# -*- coding: utf-8 -*-
"""Css styles for reports"""
import pathlib
import shutil
import logging

log = logging.getLogger(__name__)

def get_path(stylename):
    return pathlib.Path(__file__).parent/(stylename)

def copy_style(style, dest):
    p = str(get_path(style))
    log.debug("Copying style from %s to %s", p, dest)
    shutil.copy2(p, dest)
