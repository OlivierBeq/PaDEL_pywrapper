# -*- coding: utf-8

"""Python wrapper for PaDEL descriptors"""

import os
from typing import Iterable, List

import more_itertools
import pandas as pd
from bounded_pool_executor import BoundedProcessPoolExecutor
from rdkit import Chem

from . import descriptors as descriptor_types


class PaDEL:
    """PaDEL wrapper to obtain molecular descriptors."""

    def __init__(self, descriptors, ignore_3D: bool = True) -> None:
        """Instantiate a wrapper to calculate PaDEL molecular descriptors.

        :param descriptors: list of descriptors or fingerprints to be calculated
        :param ignore_3D: remove descriptors requiring 3D molecular coordinates from the provided list
        """
        # Ensure descriptors are actual PaDEL descriptors
        for descriptor in descriptors:
            if descriptor not in (descriptor_types.descriptors + descriptor_types._fingerprints) and \
                    descriptor.__class__ not in (descriptor_types.descriptors + descriptor_types._fingerprints):
                raise ValueError(f'descriptor {descriptor} is not a valid PaDEL descriptor.')

        self.descriptors: descriptor_types.Descriptor = []
        # Remove 3D descriptors if required
        self._cannot_parallelize = False
        for descriptor in descriptors:
            # Check type
            if type(descriptor).__name__ == 'ABCMeta':
                desc = descriptor()
            else:
                desc = descriptor
                self._cannot_parallelize = True
            # Do not add descriptor if ignore_3D and is 3D
            if not ignore_3D or not desc.is3D:
                self.descriptors.append(descriptor)

        self._ignore_3D = ignore_3D

    def calculate(self, mols: Iterable[Chem.Mol], show_banner: bool = True, njobs: int = 1,
                  chunksize: int = 100) -> pd.DataFrame:
        """Caclulate PaDEL descriptors.

        :param mols: RDKit molecules for which PaDEL descriptors should be calculated
        :param show_banner: If True, show notice on PaDEL descriptors usage
        :param njobs: number of concurrent processes
        :param chunksize: number of molecules to be processed by a process; ignored if njobs is 1
        :return: a pandas DataFrame containing all PaDEL descriptor values
        """
        if show_banner:
            self._show_banner()
        # Parallelize should need be
        if njobs > 1 and not self._cannot_parallelize:
            with BoundedProcessPoolExecutor(max_workers=njobs) as worker:
                futures = [worker.submit(self._calculate, list(chunk))
                           for chunk in more_itertools.batched(mols, chunksize)
                           ]
            return pd.concat([future.result()
                              for future in futures]
                             ).reset_index(drop=True)
        # Single process
        return self._calculate(list(mols))

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

    def _calculate(self, mols: List[Chem.Mol]) -> pd.DataFrame:
        """Calculate PaDEL descriptors on one process.

        :param mols: RDkit molecules for which PaDEL descriptors should be calculated.
         Only the last conformer of molecules is considered.
        :return: a pandas DataFrame containing all PaDEL desciptor values
        """
        desc_results = []
        for descriptor in self.descriptors:
            # Ensure is instantiated
            if type(descriptor).__name__ == 'ABCMeta':
                desc = descriptor()
            else:
                desc = descriptor
            desc_results.append(desc.calculate(mols))
        return pd.concat(desc_results, axis=1)

    @property
    def details(self):
        """Path to the file detailing descriptors."""
        return ("Full details about the PaDEL descriptors can be found in the file located at:\n"
                f"{os.path.join(__file__, os.pardir, 'PaDEL-Descriptor', 'Descriptors.xls')}")
