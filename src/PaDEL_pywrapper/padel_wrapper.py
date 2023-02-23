# -*- coding: utf-8

"""Python wrapper for PaDEL descriptors"""

import os
import io
import warnings
from copy import deepcopy
from typing import Iterable, List, Tuple, Union
from subprocess import PIPE, Popen

import more_itertools
import numpy as np
import pandas as pd
from bounded_pool_executor import BoundedProcessPoolExecutor
from rdkit import Chem
from rdkit.Chem import AllChem

from . import descriptor as descriptor_types
from .descriptor import Descriptor, Fingerprint
from .utils import install_java, mktempfile, needsHs


class PaDEL:
    """PaDEL wrapper to obtain molecular descriptors."""

    def __init__(self, descriptors: List[Union[Descriptor, Fingerprint]], ignore_3D: bool = True) -> None:
        """Instantiate a wrapper to calculate PaDEL molecular descriptors.

        :param descriptors: list of descriptors or fingerprints to be calculated
        :param ignore_3D: remove descriptors requiring 3D molecular coordinates from the provided list
        """
        # Ensure descriptors are actual PaDEL descriptors
        for descriptor in descriptors:
            if descriptor not in (descriptor_types.descriptors + descriptor_types._fingerprints) and \
                    descriptor.__class__ not in (descriptor_types.descriptors + descriptor_types._fingerprints):
                raise ValueError(f'descriptor {descriptor} is not a valid PaDEL descriptor.')

        self.descriptors = []
        self.fingerprints = []
        self._ignore_3D = ignore_3D
        self.has_3D_descriptors = False
        # Remove 3D descriptors if required
        for descriptor in descriptors:
            if isinstance(descriptor, Descriptor):
                # Do not add descriptor if ignore_3D and is 3D
                if not (ignore_3D and descriptor.is_3D):
                    self.descriptors.append(descriptor)
                if not ignore_3D and descriptor.is_3D:
                    self.has_3D_descriptors = True
            else:
                self.fingerprints.append(descriptor)
        # Keep registry of desc names to keep
        self.desc_kept = []
        for desc in self.descriptors:
            # Add all subdescs names
            self.desc_kept.extend(desc.subcomponents)

    def calculate(self, mols: Iterable[Chem.Mol], show_banner: bool = True, njobs: int = 1,
                  chunksize: int = 100) -> pd.DataFrame:
        """Calculate PaDEL descriptors.

        :param mols: RDKit molecules for which PaDEL descriptors should be calculated
        :param show_banner: If True, show notice on PaDEL descriptors usage
        :param njobs: number of concurrent processes
        :param chunksize: number of molecules to be processed by a process; ignored if njobs is 1
        :return: a pandas DataFrame containing all PaDEL descriptor values
        """
        if show_banner:
            self._show_banner()
        # Parallelize should need be
        if njobs > 1:
            with BoundedProcessPoolExecutor(max_workers=njobs) as worker:
                futures = [worker.submit(self._multiproc_calculate, list(chunk))
                           for chunk in more_itertools.batched(mols, chunksize)
                           ]
            # Collect results
            result = pd.concat([future.result() for future in futures]).reset_index(drop=True)
        else:
            # Single process
            result = self._calculate(list(mols))
        return result

    def _show_banner(self):
        """Print info message for citing."""
        print("""PaDEL-Descriptor is a software for calculating molecular
descriptors and fingerprints. The software calculates
1875 descriptors (1444 1D and 2D descriptors, and 431
3D descriptors) and 12 types of fingerprints.

###################################

Should you publish results based on the PaDEL descriptors,
please cite:

Yap, C.W. (2011), PaDEL-descriptor: An open source software
to calculate molecular descriptors and fingerprints.
J. Comput. Chem., 32: 1466-1474. https://doi.org/10.1002/jcc.21707

###################################
""")

    def _prepare_command(self, mols: List[Chem.Mol]) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
        """Create the ePaDEL command to be run to obtain molecular descriptors.

        :param mols: molecules to obtained molecular descriptors of
        :return: The command to run.
        """
        # 1) Ensure JRE is accessible
        self._java_path = install_java()
        # 2) Create temp SD v2k file
        self._tmp_sd = mktempfile('molecules_v2k.sd')
        self._skipped = []
        try:
            writer = Chem.SDWriter(self._tmp_sd)
            # Ensure V2000 as CDK cannot properly process v3000
            writer.SetForceV3000(False)
            for i, mol in enumerate(mols):
                if mol is not None and isinstance(mol, Chem.Mol):
                    if mol.GetNumAtoms() > 999:
                        raise ValueError('Cannot calculate descriptors for molecules with more than 999 atoms.')
                    # Does molecule lack hydrogen atoms?
                    if needsHs(mol):
                        warnings.warn('Molecule lacks hydrogen atoms: this will affect the value of calculated descriptors')
                    # If molecule has no conformer
                    confs = list(mol.GetConformers())
                    if not (len(confs) > 0 and confs[-1].Is3D()):
                        if self.has_3D_descriptors:
                            raise ValueError('Cannot calculate descriptors for a conformer-less molecule')
                        # If no 3D descriptor, compute 2D coords
                        AllChem.Compute2DCoords(mol)
                    writer.write(mol)
                else:
                    self._skipped.append(i)
            writer.close()
        except ValueError as e:
            # Free resources and raise error
            writer.close()
            os.remove(self._tmp_sd)
            raise e from None
        # 3) Create commands
        desc_commands = []
        fp_commands = []
        java_path = install_java()
        epadel_path = os.path.abspath(os.path.join(__file__, os.pardir, 'PaDEL-Descriptor', 'lib', 'ePaDEL.jar'))
        command_prefix = f"{java_path} -Djava.awt.headless=true -jar {epadel_path}"
        # Create commands for descriptors
        if len(self.descriptors):
            command_names = f"{command_prefix} -d --names"
            command_values = f"{command_prefix} -d -i {self._tmp_sd}"
            if not self._ignore_3D:
                command_names += ' -3D'
                command_values += ' -3D'
            desc_commands.append((command_names, command_values))
        # Create commands for fingerprints
        if len(self.fingerprints):
            for fp in self.fingerprints:
                command_names = f"{command_prefix} -f {fp.bit_prefix} --names"
                command_values = f"{command_prefix} -f {fp.bit_prefix} -i {self._tmp_sd}"
                # Add additional parameters
                if hasattr(fp, 'searchDepth'):
                    command_names += f" -nBits {fp.nBits} -searchDepth {fp.searchDepth}"
                    command_values += f" -nBits {fp.nBits} -searchDepth {fp.searchDepth}"
                fp_commands.append((command_names, command_values))
        return desc_commands, fp_commands

    def _cleanup(self) -> None:
        """Cleanup resources used for calculation."""
        # Remove temporary file
        os.remove(self._tmp_sd)

    def _run_command(self, commands: Tuple[str, str]) -> pd.DataFrame:
        """Run the ePaDEL command couple.

        :param commands: couple of command to be run (names and values).
        """
        with Popen(commands[0].split(), stdout=PIPE) as process:
            names = process.stdout.read().decode().split()
        with Popen(commands[1].split(), stdout=PIPE) as process:
            values = io.StringIO(process.stdout.read().decode())
            values = pd.read_csv(values,
                                 sep=' ', header=None, names=names)
        return values



    def _calculate(self, mols: List[Chem.Mol]) -> pd.DataFrame:
        """Calculate PaDEL descriptors on one process.

        :param mols: RDkit molecules for which PaDEL descriptors should be calculated.
         Only the last conformer of molecules is considered.
        :return: a pandas DataFrame containing all PaDEL desciptor values and the path to the temp dir to be removed
        """
        results = []
        # Prepare inputs
        desc_commands, fp_commands = self._prepare_command(mols)
        # If descriptors to be calculated
        if len(desc_commands):
            desc_values = self._run_command(desc_commands[0])
            # Keep only specified descriptors
            selected_features = desc_values.columns[~desc_values.columns.isin(self.desc_kept)]
            desc_values.drop(columns=selected_features, inplace=True)
            results.append(desc_values)
        # If FPs to be calculated
        for fp_command in fp_commands:
            fp_values = self._run_command(fp_command)
            results.append(fp_values)
        # Cleanup
        self._cleanup()
        results = pd.concat(results, axis=1)
        # Insert lines of skipped molecules
        if len(self._skipped):
            results = pd.DataFrame(np.insert(results.values, self._skipped,
                                             values=[np.NaN] * len(results.columns),
                                             axis=0),
                                   columns=results.columns)
        return results


    def _multiproc_calculate(self, mols: List[Chem.Mol]) -> pd.DataFrame:
        """Calculate PaDEL descriptors in thread-safe manner.

        :param mols: RDkit molecules for which PaDEL descriptors should be calculated.
         Only the last conformer of molecules is considered.
        :return: a pandas DataFrame containing all PaDEL desciptor values and the path to the temp dir to be removed
        """
        # Copy self instance to make thread safe
        padel = deepcopy(self)
        # Run copy
        result = padel.calculate(mols, show_banner=False, njobs=1)
        return result


    @property
    def details(self):
        """Path to the file detailing descriptors."""
        return ("Full details about the PaDEL descriptors can be found in the file located at:\n"
                f"{os.path.abspath(os.path.join(__file__, os.pardir, 'PaDEL-Descriptor', 'Descriptors.xls'))}")
