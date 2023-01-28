# -*- coding: utf-8

"""Utility functions."""

from typing import Union

import numpy as np


def parse_numeric(value: str) -> Union[int, float]:
    """Parse a string representation of a number."""
    if not len(value):
        return np.nan
    return int(value) if value.lstrip(" -+").isdigit() else float(value)


def parse_numeric_to_null(value: str) -> Union[int, float]:
    """Parse a string representation of a number and return NaN."""
    return np.nan
