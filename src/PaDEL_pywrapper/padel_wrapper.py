"""Python wrapper for PaDEL descriptors"""

import io
import multiprocessing
import os
import subprocess
import warnings
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy

import more_itertools
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from . import descriptor as descriptor_types
from .descriptor import Descriptor, Fingerprint
from .utils import install_java, mktempfile, needsHs


def _make_chunks(mols: list, njobs: int, chunksize: int | None) -> list[list]:
    """Split molecules into chunks to be dispatched to worker processes.

    If chunksize is None, mols are split into exactly min(njobs, len(mols))
    balanced chunks (sizes differ by at most one) so every requested worker
    is guaranteed to receive work. If chunksize is given, it is authoritative
    and mols are batched into fixed-size chunks (may leave some workers idle).

    :param mols: molecules (or any list) to split into chunks
    :param njobs: number of concurrent workers that will consume the chunks
    :param chunksize: fixed size of each chunk, or None to auto-balance across njobs
    :return: a list of chunks
    """
    if not mols:
        return []
    if chunksize is None:
        n_chunks = min(njobs, len(mols))
        base, remainder = divmod(len(mols), n_chunks)
        chunks = []
        start = 0
        for i in range(n_chunks):
            size = base + (1 if i < remainder else 0)
            chunks.append(mols[start:start + size])
            start += size
        return chunks
    return [list(chunk) for chunk in more_itertools.batched(mols, chunksize)]


class PaDEL:
    """PaDEL wrapper to obtain molecular descriptors."""

    # Extra safety net for JRE installation; the real guarantee is that calculate() installs
    # the JRE once, before any worker starts.
    lock = multiprocessing.RLock()

    def __init__(self, descriptors: list[Descriptor | Fingerprint], ignore_3D: bool = True) -> None:
        """Instantiate a wrapper to calculate PaDEL molecular descriptors.

        :param descriptors: list of descriptors or fingerprints to be calculated
        :param ignore_3D: remove descriptors requiring 3D molecular coordinates from the provided list
        """
        # Ensure descriptors are actual PaDEL descriptors
        names = ([descriptor.name for descriptor in descriptor_types.descriptors] +
                 [fingerprint.name for fingerprint in descriptor_types._fingerprints])
        for descriptor in descriptors:
            if not isinstance(descriptor, (Descriptor, Fingerprint)) or descriptor.name not in names:
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
                  chunksize: int | None = None, timeout: float | None = None) -> pd.DataFrame:
        """Calculate PaDEL descriptors.

        :param mols: RDKit molecules for which PaDEL descriptors should be calculated
        :param show_banner: If True, show notice on PaDEL descriptors usage
        :param njobs: number of concurrent processes; must not exceed the number of available
            CPU cores. Parallelism comes from spawning njobs OS processes, each running its own
            single-core-pinned JVM (-XX:ActiveProcessorCount=1), so JVMs do not oversubscribe
            the host by each auto-sizing GC/JIT thread pools to the full core count.
        :param chunksize: number of molecules processed per worker process; if None (default)
            molecules are auto-balanced across njobs workers so every worker gets work; ignored
            if njobs is 1
        :param timeout: maximum number of seconds allowed for each descriptor/fingerprint
            calculation for a single molecule. None (default) waits indefinitely. This is
            enforced on the Java side (see ``-timeout``): a descriptor exceeding it for one
            molecule returns NaN for that molecule instead of hanging the whole batch. The
            underlying ePaDEL subprocess call itself is also bounded, by a multiple of this
            value large enough to cover every (molecule, descriptor/fingerprint) pair it could
            legitimately need to wait out, so that outer safety net cannot fire before Java's
            own per-item timeout has had a chance to resolve all of them. Some PaDEL 3D
            descriptors (e.g. WHIM) rely on an eigenvalue decomposition that can fail to
            converge and hang indefinitely for molecules with degenerate 3D geometry (very few
            heavy atoms, near-linear/near-planar arrangements); setting a timeout guards
            against this instead of hanging forever.
        :return: a pandas DataFrame containing all PaDEL descriptor values
        """
        if njobs < 1:
            raise ValueError('njobs must be a strictly positive integer.')
        cpu_count = os.cpu_count() or 1
        if njobs > cpu_count:
            raise ValueError(f'njobs ({njobs}) exceeds the number of available CPU cores ({cpu_count}).')
        if show_banner:
            self._show_banner()
        mols = list(mols)
        # Parallelize should need be
        if njobs > 1:
            # Ensure the JRE is installed once in the parent process before workers start,
            # so concurrently-spawned workers all find it already in place.
            with self.lock:
                install_java()
            chunks = _make_chunks(mols, njobs, chunksize)
            if not chunks:
                return pd.DataFrame()
            with ProcessPoolExecutor(max_workers=njobs) as worker:
                futures = [worker.submit(self._multiproc_calculate, chunk, timeout) for chunk in chunks]
                results = [future.result() for future in futures]
            # Collect results
            result = pd.concat(results).reset_index(drop=True)
        else:
            # Single process
            result = self._calculate(mols, timeout=timeout)
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

    def _prepare_command(self, mols: list[Chem.Mol], timeout: float | None = None
                          ) -> tuple[list[tuple[list[str], list[str]]],
                                     list[tuple[list[str], list[str]]]]:
        """Create the ePaDEL command to be run to obtain molecular descriptors.

        :param mols: molecules to obtained molecular descriptors of
        :param timeout: maximum number of seconds ePaDEL is allowed to spend calculating a
            single descriptor or fingerprint for a single molecule; forwarded to the Java
            process as ``-timeout``, so that e.g. WHIM's eigenvalue decomposition failing to
            converge on one molecule returns NaN for it instead of hanging the whole batch.
            None (default) does not pass ``-timeout``, so Java waits indefinitely.
        :return: The commands (as argument lists) to run.
        """
        # 1) Ensure JRE is accessible
        with self.lock:
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
                        warnings.warn('Molecule lacks hydrogen atoms: this will affect the value '
                                      'of calculated descriptors', stacklevel=2)
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
        epadel_path = os.path.abspath(os.path.join(__file__, os.pardir, 'PaDEL-Descriptor', 'lib', 'ePaDEL.jar'))
        # Pin each JVM to a single core: parallelism comes from spawning njobs OS
        # processes, not from letting each JVM auto-size GC/JIT threads to the host's
        # full core count (which would oversubscribe cores beyond what njobs requested).
        command_prefix = [self._java_path, '-Djava.awt.headless=true',
                          '-XX:ActiveProcessorCount=1', '-jar', epadel_path]
        # Create commands for descriptors
        if len(self.descriptors):
            command_names = [*command_prefix, '-d', '--names']
            command_values = [*command_prefix, '-d', '-i', self._tmp_sd]
            if not self._ignore_3D:
                command_names = [*command_names, '-3D']
                command_values = [*command_values, '-3D']
            if timeout is not None:
                command_values = [*command_values, '-timeout', str(timeout)]
            desc_commands.append((command_names, command_values))
        # Create commands for fingerprints
        if len(self.fingerprints):
            for fp in self.fingerprints:
                command_names = [*command_prefix, '-f', fp.bit_prefix, '--names']
                command_values = [*command_prefix, '-f', fp.bit_prefix, '-i', self._tmp_sd]
                # Add additional parameters
                if fp.size is not None:
                    extra = ['-nBits', str(fp.nBits), '-searchDepth', str(fp.searchDepth)]
                    command_names = [*command_names, *extra]
                    command_values = [*command_values, *extra]
                if timeout is not None:
                    command_values = [*command_values, '-timeout', str(timeout)]
                fp_commands.append((command_names, command_values))
        return desc_commands, fp_commands

    def _cleanup(self) -> None:
        """Cleanup resources used for calculation."""
        # Remove temporary file
        os.remove(self._tmp_sd)

    def _run_command(self, commands: tuple[list[str], list[str]], expected_rows: int,
                      timeout: float | None = None) -> pd.DataFrame:
        """Run the ePaDEL command couple.

        :param commands: A couple of commands to be run (names and values).
        :param expected_rows: number of molecules actually written to the SD file
            (len(mols) - len(self._skipped)); used to detect a truncated/crashed
            Java process instead of silently returning misaligned data.
        :param timeout: maximum number of seconds to wait for each subprocess; None waits
            indefinitely. See :meth:`calculate` for why this matters (PCA-based 3D
            descriptors like WHIM can hang on degenerate geometries).
        :raises RuntimeError: if either subprocess exits non-zero, times out, or if the
            number of value rows returned does not match expected_rows.
        """
        try:
            names_proc = subprocess.run(commands[0], capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired as e:
            raise RuntimeError(
                f"ePaDEL timed out after {timeout}s while fetching descriptor/fingerprint names."
            ) from e
        if names_proc.returncode != 0:
            raise RuntimeError(
                f"ePaDEL failed while fetching descriptor/fingerprint names "
                f"(exit code {names_proc.returncode}):\n{names_proc.stderr.decode(errors='replace')}"
            )
        names = names_proc.stdout.decode().split()

        try:
            values_proc = subprocess.run(commands[1], capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired as e:
            raise RuntimeError(
                f"ePaDEL timed out after {timeout}s while calculating values. This can happen "
                f"for PCA-based 3D descriptors (e.g. WHIM) on molecules with degenerate 3D "
                f"geometry (very few heavy atoms, near-linear/near-planar arrangements), whose "
                f"eigenvalue decomposition can fail to converge."
            ) from e
        if values_proc.returncode != 0:
            raise RuntimeError(
                f"ePaDEL failed while calculating values "
                f"(exit code {values_proc.returncode}):\n{values_proc.stderr.decode(errors='replace')}"
            )
        values = pd.read_csv(io.StringIO(values_proc.stdout.decode()), sep=' ', header=None, names=names)

        if len(values) != expected_rows:
            raise RuntimeError(
                f"ePaDEL returned {len(values)} rows but {expected_rows} were expected "
                f"(command likely crashed partway through, silently truncating output). "
                f"stderr:\n{values_proc.stderr.decode(errors='replace')}"
            )
        return values

    @staticmethod
    def _subprocess_timeout(per_item_timeout: float | None, n_mols: int, n_items: int) -> float | None:
        """Derive the outer ePaDEL subprocess timeout from the per-descriptor/fingerprint
        Java-side timeout.

        Java bounds each individual (molecule, descriptor/fingerprint) calculation to
        ``per_item_timeout`` seconds (see ``-timeout`` in :meth:`_prepare_command`), so the
        whole subprocess call may legitimately take up to ``n_mols * n_items`` times that long
        (e.g. one molecule with degenerate 3D geometry making a single 3D descriptor time out).
        The outer subprocess-level timeout must therefore be scaled accordingly -- passing
        ``per_item_timeout`` straight through would let the outer timeout fire before Java's own
        per-item timeout had a chance to resolve every pair, defeating the point of the
        per-item mechanism.

        :param per_item_timeout: seconds allotted to each (molecule, descriptor/fingerprint)
            pair on the Java side; None means no timeout at all.
        :param n_mols: number of molecules in this call.
        :param n_items: number of descriptors or fingerprint types computed by this call.
        :return: the outer subprocess-level timeout in seconds, or None.
        """
        if per_item_timeout is None:
            return None
        return per_item_timeout * max(n_mols, 1) * max(n_items, 1)

    def _calculate(self, mols: list[Chem.Mol], timeout: float | None = None) -> pd.DataFrame:
        """Calculate PaDEL descriptors on one process.

        :param mols: RDkit molecules for which PaDEL descriptors should be calculated.
         Only the last conformer of molecules is considered.
        :param timeout: maximum number of seconds allowed for each descriptor/fingerprint
            calculation for a single molecule, on the Java side; None waits indefinitely.
            See :meth:`calculate`.
        :return: a pandas DataFrame containing all PaDEL desciptor values and the path to the temp dir to be removed
        """
        results = []
        # Prepare inputs
        desc_commands, fp_commands = self._prepare_command(mols, timeout)
        expected_rows = len(mols) - len(self._skipped)
        # If descriptors to be calculated
        if len(desc_commands):
            # Java always computes the full descriptor catalog (Python only subsets the
            # columns afterward), so the outer timeout must account for every one of them.
            n_desc_items = sum(1 for d in descriptor_types.descriptors if not self._ignore_3D or not d.is_3D)
            desc_values = self._run_command(desc_commands[0], expected_rows,
                                             self._subprocess_timeout(timeout, len(mols), n_desc_items))
            # Keep only specified descriptors
            selected_features = desc_values.columns[~desc_values.columns.isin(self.desc_kept)]
            desc_values.drop(columns=selected_features, inplace=True)
            results.append(desc_values)
        # If FPs to be calculated
        for fp_command in fp_commands:
            # Each command computes a single fingerprint type across all molecules.
            fp_values = self._run_command(fp_command, expected_rows,
                                           self._subprocess_timeout(timeout, len(mols), 1))
            results.append(fp_values)
        # Cleanup
        self._cleanup()
        results = pd.concat(results, axis=1)
        # Insert lines of skipped molecules
        if len(self._skipped):
            results = pd.DataFrame(np.insert(results.values, self._skipped,
                                             values=[np.nan] * len(results.columns),
                                             axis=0),
                                   columns=results.columns)
        return results

    def _multiproc_calculate(self, mols: list[Chem.Mol], timeout: float | None = None) -> pd.DataFrame:
        """Calculate PaDEL descriptors in a worker process.

        :param mols: RDkit molecules for which PaDEL descriptors should be calculated.
         Only the last conformer of molecules is considered.
        :param timeout: maximum number of seconds to wait for each underlying ePaDEL process;
            None waits indefinitely. See :meth:`calculate`.
        :return: a pandas DataFrame containing all PaDEL desciptor values and the path to the temp dir to be removed
        """
        # Copy self instance to make thread safe
        padel = deepcopy(self)
        # Run copy
        result = padel.calculate(mols, show_banner=False, njobs=1, timeout=timeout)
        return result

    @property
    def details(self):
        """Path to the file detailing descriptors."""
        return ("Full details about the PaDEL descriptors can be found in the file located at:\n"
                f"{os.path.abspath(os.path.join(__file__, os.pardir, 'PaDEL-Descriptor', 'Descriptors.xls'))}")
