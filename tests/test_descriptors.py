"""Tests for :class:`PaDEL_pywrapper.descriptor.Descriptor` and molecular descriptor calculation."""

import contextlib
import io
import unittest

from PaDEL_pywrapper import PaDEL
from PaDEL_pywrapper.descriptor import ALOGP, TPSA, WHIM, Descriptor, Fingerprint, GraphOnlyFP, Weight
from tests.constants import MOLECULES


class TestDescriptorClass(unittest.TestCase):
    """Tests for the Descriptor class metadata loaded from padel_descs.tsv."""

    def test_descriptor_subcomponents_populated(self):
        self.assertIn('ALogP', ALOGP.subcomponents)
        self.assertTrue(len(ALOGP.subcomponents) > 0)

    def test_descriptor_description_populated(self):
        self.assertFalse(ALOGP.description.empty)
        self.assertListEqual(list(ALOGP.description.columns), ['name', 'description'])

    def test_descriptor_is_3d_flag(self):
        self.assertFalse(ALOGP.is_3D)
        self.assertTrue(WHIM.is_3D)

    def test_unknown_descriptor_name_has_no_subcomponents(self):
        bogus = Descriptor('TotallyBogusName', False)
        self.assertEqual(bogus.subcomponents, [])


class TestConfigurableFingerprint(unittest.TestCase):
    """Tests for Fingerprint.__call__'s size/searchDepth configuration."""

    def test_non_configurable_fingerprint_size_stays_none(self):
        fp = Fingerprint('MACCS fingerprint')
        self.assertIsNone(fp.size)
        self.assertIsNone(fp.searchDepth)

    def test_call_sets_size_and_search_depth_for_configurable_fingerprint(self):
        custom = GraphOnlyFP(size=2048, searchDepth=8)
        self.assertEqual(custom.size, 2048)
        self.assertEqual(custom.nBits, 2048)
        self.assertEqual(custom.searchDepth, 8)
        # The original singleton instance is untouched.
        self.assertNotEqual(GraphOnlyFP.size, 2048)


class TestDescriptorCalculation(unittest.TestCase):
    """Tests for calculating real 1D/2D/3D descriptors end-to-end."""

    def setUp(self) -> None:
        # Cisplatin triggers a known upstream Java ClassCastException in EStateAtomType
        # (see tests/test_error_paths.py) -- excluded here to test the happy path.
        self.molecules = [mol for name, mol in MOLECULES.items() if name != 'cisplatin']

    def test_2d_descriptor_calculation(self):
        padel = PaDEL([ALOGP, Weight, TPSA])
        values = padel.calculate(self.molecules, show_banner=False)
        self.assertEqual(values.shape[0], len(self.molecules))
        self.assertEqual(set(values.columns), set(padel.desc_kept))
        self.assertFalse(values.isna().any().any())

    def test_3d_descriptor_calculation(self):
        padel = PaDEL([WHIM], ignore_3D=False)
        self.assertTrue(padel.has_3D_descriptors)
        values = padel.calculate(self.molecules, show_banner=False)
        self.assertEqual(values.shape[0], len(self.molecules))
        self.assertFalse(values.isna().any().any())

    def test_3d_descriptor_dropped_when_ignore_3d(self):
        padel = PaDEL([WHIM], ignore_3D=True)
        self.assertEqual(padel.descriptors, [])
        self.assertFalse(padel.has_3D_descriptors)

    def test_show_banner_prints_citation_notice(self):
        padel = PaDEL([ALOGP])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            padel.calculate([self.molecules[0]], show_banner=True)
        self.assertIn('PaDEL-Descriptor', buf.getvalue())

    def test_details_property_points_to_descriptors_file(self):
        padel = PaDEL([ALOGP])
        self.assertIn('Descriptors.xls', padel.details)


if __name__ == "__main__":
    unittest.main()
