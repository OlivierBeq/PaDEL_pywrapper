"""Regression tests for PaDEL descriptor/fingerprint validation."""

import unittest

from PaDEL_pywrapper import PaDEL
from PaDEL_pywrapper.descriptor import ALOGP, FP, Descriptor


class TestValidation(unittest.TestCase):
    """Tests for PaDEL.__init__'s validation of the descriptors/fingerprints list."""

    def test_rejects_bogus_named_descriptor_instance(self):
        """A Descriptor instance with an unrecognized name must be rejected.

        Regression test: the original validation used `and` instead of `or`, so any
        actual Descriptor/Fingerprint instance bypassed the name check entirely.
        """
        bogus = Descriptor('TotallyBogusName', False)
        with self.assertRaises(ValueError):
            PaDEL([bogus])

    def test_rejects_non_descriptor_object(self):
        with self.assertRaises(ValueError):
            PaDEL(['not_a_descriptor_object'])

    def test_accepts_valid_descriptor_and_fingerprint(self):
        padel = PaDEL([ALOGP, FP])
        self.assertIn(ALOGP, padel.descriptors)
        self.assertIn(FP, padel.fingerprints)


if __name__ == "__main__":
    unittest.main()
