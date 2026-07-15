"""Tests for error paths: oversized molecules, missing conformers, skipped/None
molecules, and detection of a crashed/truncated ePaDEL subprocess."""

import subprocess
import unittest
from unittest.mock import MagicMock, patch

from rdkit import Chem

from PaDEL_pywrapper import PaDEL
from PaDEL_pywrapper.descriptor import ALOGP, WHIM, descriptors
from tests.constants import MOLECULES


class TestOversizedMolecule(unittest.TestCase):
    def test_raises_for_more_than_999_atoms(self):
        mol = MOLECULES['erlotinib']
        with patch.object(Chem.Mol, 'GetNumAtoms', return_value=1000):
            with self.assertRaises(ValueError):
                PaDEL([ALOGP]).calculate([mol], show_banner=False)


class TestConformerlessMolecule(unittest.TestCase):
    def test_raises_for_conformerless_molecule_with_3d_requested(self):
        mol = Chem.AddHs(Chem.MolFromSmiles('CCC'))  # no conformer embedded
        padel = PaDEL([WHIM], ignore_3D=False)
        with self.assertRaises(ValueError):
            padel.calculate([mol], show_banner=False)

    def test_conformerless_molecule_ok_without_3d_descriptors(self):
        mol = Chem.AddHs(Chem.MolFromSmiles('CCC'))
        padel = PaDEL([ALOGP])
        # Should not raise: 2D coordinates are computed instead.
        values = padel.calculate([mol], show_banner=False)
        self.assertEqual(values.shape[0], 1)


class TestMissingHydrogensWarning(unittest.TestCase):
    def test_warns_when_molecule_lacks_explicit_hydrogens(self):
        mol = Chem.MolFromSmiles('CCC')  # no explicit Hs added
        padel = PaDEL([ALOGP])
        with self.assertWarns(UserWarning):
            padel.calculate([mol], show_banner=False)


class TestSkippedMolecules(unittest.TestCase):
    def test_none_molecule_skipped_and_reinserted_as_nan_row(self):
        mols = [MOLECULES['erlotinib'], None, MOLECULES['midecamycin']]
        result = PaDEL([ALOGP]).calculate(mols, show_banner=False)
        self.assertEqual(result.shape[0], 3)
        self.assertTrue(result.iloc[1].isna().all())
        self.assertFalse(result.iloc[0].isna().any())
        self.assertFalse(result.iloc[2].isna().any())


class TestTruncatedSubprocessDetection(unittest.TestCase):
    """Mocked subprocess tests for the row-count/returncode safety net, independent
    of which ePaDEL.jar is installed."""

    def test_nonzero_returncode_raises_runtime_error(self):
        padel = PaDEL([ALOGP])
        mols = [MOLECULES['erlotinib']]
        names_result = MagicMock(returncode=0, stdout=b'A B C', stderr=b'')
        values_result = MagicMock(returncode=1, stdout=b'', stderr=b'boom')
        with patch('PaDEL_pywrapper.padel_wrapper.subprocess.run', side_effect=[names_result, values_result]):
            with self.assertRaises(RuntimeError):
                padel.calculate(mols, show_banner=False)

    def test_nonzero_returncode_fetching_names_raises_runtime_error(self):
        padel = PaDEL([ALOGP])
        mols = [MOLECULES['erlotinib']]
        names_result = MagicMock(returncode=1, stdout=b'', stderr=b'names boom')
        with patch('PaDEL_pywrapper.padel_wrapper.subprocess.run', return_value=names_result):
            with self.assertRaises(RuntimeError):
                padel.calculate(mols, show_banner=False)

    def test_row_count_mismatch_raises_runtime_error(self):
        padel = PaDEL([ALOGP])
        mols = [MOLECULES['erlotinib'], MOLECULES['midecamycin']]  # expected_rows == 2
        names_result = MagicMock(returncode=0, stdout=b'A B C', stderr=b'')
        # Only one row of values instead of the expected two: simulates a JVM
        # crashing partway through the batch.
        values_result = MagicMock(returncode=0, stdout=b'1 2 3\n', stderr=b'crashed mid-batch')
        with patch('PaDEL_pywrapper.padel_wrapper.subprocess.run', side_effect=[names_result, values_result]):
            with self.assertRaises(RuntimeError):
                padel.calculate(mols, show_banner=False)

    def test_timeout_fetching_names_raises_runtime_error(self):
        padel = PaDEL([ALOGP])
        mols = [MOLECULES['erlotinib']]
        with patch('PaDEL_pywrapper.padel_wrapper.subprocess.run',
                  side_effect=subprocess.TimeoutExpired(cmd='java', timeout=5)):
            with self.assertRaises(RuntimeError):
                padel.calculate(mols, show_banner=False, timeout=5)

    def test_timeout_calculating_values_raises_runtime_error(self):
        # Simulates a hung WHIM-style eigenvalue decomposition: names succeed, values hang.
        padel = PaDEL([ALOGP])
        mols = [MOLECULES['erlotinib']]
        names_result = MagicMock(returncode=0, stdout=b'A B C', stderr=b'')
        with patch('PaDEL_pywrapper.padel_wrapper.subprocess.run',
                  side_effect=[names_result, subprocess.TimeoutExpired(cmd='java', timeout=5)]):
            with self.assertRaises(RuntimeError):
                padel.calculate(mols, show_banner=False, timeout=5)


class TestCisplatinNeverSilentlyMisaligned(unittest.TestCase):
    """Documents the current behavior for the known Pt/cisplatin ClassCastException in
    EStateAtomType, plus a second, previously-masked failure mode it exposed once the
    Main.java catch was broadened: WHIM's PCA eigenvalue decomposition (Jama's `hqr2`
    Hessenberg QR algorithm) does not converge for cisplatin's degenerate 3D geometry and
    hangs indefinitely instead of raising. With the old jar this was never reached because
    EStateAtomType's uncaught ClassCastException killed the process first; with the
    broadened catch, execution now reaches WHIM and hangs. Whichever jar is installed, the
    invariants "never silently misaligned" and "never hangs forever" must both hold -- the
    latter is enforced here via the `timeout` safety net."""

    def test_cisplatin_full_descriptor_set_never_hangs_or_misaligns(self):
        mols = [MOLECULES['erlotinib'], MOLECULES['cisplatin']]
        padel = PaDEL(descriptors, ignore_3D=False)
        try:
            df = padel.calculate(mols, show_banner=False, timeout=20)
        except RuntimeError:
            # Either the JVM crashed on cisplatin's Pt atom (row-count/returncode safety
            # net) or WHIM's eigenvalue decomposition failed to converge in time (timeout
            # safety net). Either way: raised loudly, never a silent hang or misalignment.
            return
        # If the jar's WHIM/Jama implementation is ever fixed upstream: all rows present,
        # with NaNs only in the specific failing sub-descriptor columns.
        self.assertEqual(df.shape[0], len(mols))


if __name__ == "__main__":
    unittest.main()
