# -*- coding: utf-8 -*-
"""Tests for molecular fingerprints."""

import unittest

from PaDEL_pywrapper import PaDEL, descriptors
from PaDEL_pywrapper.descriptor import _fingerprints
from tests.constants import MOLECULES


class TestFingerprints(unittest.TestCase):
    """Tests for PaDEL_pywrapper molecular fingerprints."""
    def setUp(self) -> None:
        """Load molecules."""
        self.molecules = list(MOLECULES.values())
        self.fp_lens = {'EStateFP': 79,
                        'MACCSFP': 166,
                        'PubchemFP': 881,
                        'SubFP': 307,
                        'SubFPC': 307,
                        'KRFP': 4860,
                        'KRFPC': 4860,
                        'AP2DFP': 780,
                        'AP2DFPC': 780
                        }

    def test_fingerprints(self):
        """Test the dimensions of the output fingerprint dataframe."""
        for fp_type in _fingerprints:
            padel = PaDEL([fp_type])
            values = padel.calculate(self.molecules, show_banner=False)
            self.assertEqual(values.shape, (len(MOLECULES), self.fp_lens.get(fp_type.short_name, 1024)))
            self.assertEqual(len(values.columns.unique().tolist()), self.fp_lens.get(fp_type.short_name, 1024))
            self.assertFalse(values.isna().any().any())

    def test_fingerprint_multithread(self):
        """Test the dimensions of the output fingerprint dataframes calculated by different processes."""
        for fp_type in _fingerprints:
            padel = PaDEL([fp_type])
            values = padel.calculate(self.molecules, show_banner=False,
                                     njobs=1, chunksize=1)
            self.assertEqual(values.shape, (len(MOLECULES), self.fp_lens.get(fp_type.short_name, 1024)))
            self.assertEqual(len(values.columns.unique().tolist()), self.fp_lens.get(fp_type.short_name, 1024))
            self.assertFalse(values.isna().any().any())
