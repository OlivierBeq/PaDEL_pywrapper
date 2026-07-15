"""Tests for parallel dispatch: chunk sizing, core-count validation, JVM pinning."""

import unittest
from unittest.mock import patch

import pandas as pd

from PaDEL_pywrapper import PaDEL
from PaDEL_pywrapper.descriptor import MACCSFP
from PaDEL_pywrapper.padel_wrapper import _make_chunks
from tests.constants import MOLECULES


class TestMakeChunks(unittest.TestCase):
    """Unit tests for the pure chunk-splitting helper (no Java involved)."""

    def test_no_molecules(self):
        self.assertEqual(_make_chunks([], njobs=4, chunksize=None), [])

    def test_auto_split_uses_all_workers_when_enough_molecules(self):
        mols = list(range(7))
        chunks = _make_chunks(mols, njobs=4, chunksize=None)
        self.assertEqual(len(chunks), 4)
        self.assertEqual(sorted(len(c) for c in chunks), [1, 2, 2, 2])
        self.assertEqual(sum(chunks, []), mols)

    def test_auto_split_fewer_molecules_than_workers(self):
        mols = list(range(3))
        chunks = _make_chunks(mols, njobs=8, chunksize=None)
        # Cannot use more workers than there are molecules: exactly len(mols) chunks.
        self.assertEqual(len(chunks), 3)
        self.assertEqual(sorted(len(c) for c in chunks), [1, 1, 1])

    def test_explicit_chunksize_is_authoritative(self):
        mols = list(range(7))
        chunks = _make_chunks(mols, njobs=4, chunksize=1000)
        # A too-large explicit chunksize still yields a single chunk (old batching behavior).
        self.assertEqual(len(chunks), 1)
        self.assertEqual(chunks[0], mols)

    def test_explicit_small_chunksize(self):
        mols = list(range(7))
        chunks = _make_chunks(mols, njobs=4, chunksize=3)
        self.assertEqual([len(c) for c in chunks], [3, 3, 1])


class TestNjobsValidation(unittest.TestCase):
    """Tests for njobs bounds checking, without invoking Java."""

    def test_njobs_zero_raises(self):
        padel = PaDEL([MACCSFP])
        with self.assertRaises(ValueError):
            padel.calculate(list(MOLECULES.values())[:2], njobs=0, show_banner=False)

    def test_njobs_exceeds_cpu_count_raises(self):
        padel = PaDEL([MACCSFP])
        with patch("os.cpu_count", return_value=2):
            with self.assertRaises(ValueError):
                padel.calculate(list(MOLECULES.values())[:2], njobs=3, show_banner=False)


class TestJvmCommandPinning(unittest.TestCase):
    """Tests that every Java invocation is pinned to a single core."""

    def test_prepare_command_pins_processor_count(self):
        padel = PaDEL([MACCSFP])
        mols = list(MOLECULES.values())[:2]
        desc_commands, fp_commands = padel._prepare_command(mols)
        self.addCleanup(padel._cleanup)
        for names_argv, values_argv in [*desc_commands, *fp_commands]:
            self.assertIn('-XX:ActiveProcessorCount=1', names_argv)
            self.assertIn('-XX:ActiveProcessorCount=1', values_argv)


class TestParallelMatchesSingleProcess(unittest.TestCase):
    """Ensures the njobs>1 code path is exercised and produces the same results as njobs=1."""

    def test_parallel_matches_single_process(self):
        mols = [mol for name, mol in MOLECULES.items() if name != 'cisplatin']
        single = PaDEL([MACCSFP]).calculate(mols, show_banner=False, njobs=1)
        parallel = PaDEL([MACCSFP]).calculate(mols, show_banner=False, njobs=2, chunksize=None)
        pd.testing.assert_frame_equal(single, parallel)

    def test_parallel_uses_all_requested_workers(self):
        # With 6 molecules and njobs=3, _make_chunks must produce 3 chunks: verified directly
        # (avoids depending on OS-level process introspection, which is flaky in CI).
        mols = [mol for name, mol in MOLECULES.items() if name != 'cisplatin']
        chunks = _make_chunks(mols, njobs=3, chunksize=None)
        self.assertEqual(len(chunks), 3)

    def test_njobs_greater_than_one_with_no_molecules_returns_empty_dataframe(self):
        result = PaDEL([MACCSFP]).calculate([], show_banner=False, njobs=2)
        self.assertTrue(result.empty)

    def test_multiproc_calculate_matches_direct_calculate(self):
        # Exercise _multiproc_calculate directly (in-process) rather than through the
        # ProcessPoolExecutor, so its body is measured by coverage.
        mols = [mol for name, mol in MOLECULES.items() if name != 'cisplatin'][:2]
        padel = PaDEL([MACCSFP])
        via_multiproc = padel._multiproc_calculate(mols)
        direct = PaDEL([MACCSFP]).calculate(mols, show_banner=False)
        pd.testing.assert_frame_equal(via_multiproc, direct)


if __name__ == "__main__":
    unittest.main()
