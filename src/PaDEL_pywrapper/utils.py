# -*- coding: utf-8

"""Utility functions."""

import sys
import os
from typing import Union
import glob

import numpy as np
from jdk import install as _jre_install, _JRE_DIR


def parse_numeric(value: str) -> Union[int, float]:
    """Parse a string representation of a number."""
    if not len(value):
        return np.nan
    return int(value) if value.lstrip(" -+").isdigit() else float(value)


def parse_numeric_to_null(value: str) -> Union[int, float]:
    """Parse a string representation of a number and return NaN."""
    return np.nan


def install_jre(version: int = 11):
    """Install a Java Runtime Environment."""
    path = glob.glob(os.path.join(_JRE_DIR, '**', 'server',
                                  'jvm.dll' if sys.platform == "win32" else 'libjvm.so'
                                 ), recursive=True)
    if len(path) == 0:
        # Could not find JRE, install it
        path = _jre_install(version, jre=True)
        path = glob.glob(os.path.join(_JRE_DIR, '**', 'server',
                                      'jvm.dll' if sys.platform == "win32" else 'libjvm.so'
                                      ), recursive=True)
    path = path[0]
    return path
