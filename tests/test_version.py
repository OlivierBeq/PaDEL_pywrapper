# -*- coding: utf-8 -*-
"""Trivial version test."""

import unittest

import PaDEL_pywrapper


class TestVersion(unittest.TestCase):
    """Trivially test a version."""
    def test_version_type(self):
        """Test the version is a string.

        This is only meant to be an example test.
        """
        version = PaDEL_pywrapper.__version__
        self.assertIsInstance(version, str)
