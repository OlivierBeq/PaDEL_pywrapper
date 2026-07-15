"""Wrapper for PaDEL descriptors"""

from .descriptor import descriptors
from .padel_wrapper import PaDEL

__all__ = ["PaDEL", "descriptors"]

__version__ = "1.1.0"
