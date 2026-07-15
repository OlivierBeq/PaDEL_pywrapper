"""Tests for :mod:`PaDEL_pywrapper.utils`."""

import os
import unittest
from unittest.mock import patch

from rdkit import Chem

from PaDEL_pywrapper import utils
from PaDEL_pywrapper.utils import get_java_in_dir, install_java, mktempdir, mktempfile, needsHs


class TestNeedsHs(unittest.TestCase):
    def test_needs_hs_without_explicit_hydrogens(self):
        mol = Chem.MolFromSmiles("CCO")
        self.assertTrue(needsHs(mol))

    def test_does_not_need_hs_with_explicit_hydrogens(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        self.assertFalse(needsHs(mol))


class TestTempFileHelpers(unittest.TestCase):
    def test_mktempdir_creates_writeable_dir(self):
        path = mktempdir()
        self.addCleanup(os.rmdir, path)
        self.assertTrue(os.path.isdir(path))

    def test_mktempfile_creates_writeable_file(self):
        path = mktempfile(suffix="foo.txt")
        self.addCleanup(os.remove, path)
        self.assertTrue(os.path.isfile(path))
        self.assertTrue(path.endswith("foo.txt"))


class TestJavaInstall(unittest.TestCase):
    def test_install_java_returns_existing_path(self):
        path = install_java()
        self.assertIsNotNone(path)
        self.assertTrue(os.path.isfile(path))

    def test_get_java_in_dir_missing_version_returns_none(self):
        search_dir = mktempdir()
        self.addCleanup(os.rmdir, search_dir)
        path = get_java_in_dir(search_dir, version=999)
        self.assertIsNone(path)

    def test_install_java_triggers_download_when_not_found(self):
        # Simulate a fresh environment: get_java_in_dir first returns None (not found),
        # then a fake path once _jre_install has "installed" it.
        with patch.object(utils, 'get_java_in_dir', side_effect=[None, '/fake/jdk-11/bin/java']) as mock_get, \
             patch.object(utils, '_jre_install') as mock_install:
            path = install_java(version=11)
        self.assertEqual(path, '/fake/jdk-11/bin/java')
        mock_install.assert_called_once_with(11, jre=True)
        self.assertEqual(mock_get.call_count, 2)


if __name__ == "__main__":
    unittest.main()
