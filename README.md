# 🧪 PaDEL Python Wrapper

<!-- Badges -->
<div align="center">

[![PyPI version](https://img.shields.io/pypi/v/PaDEL-pywrapper.svg)](https://pypi.org/project/PaDEL-pywrapper/)
[![Supported Python versions](https://img.shields.io/pypi/pyversions/PaDEL-pywrapper.svg)](https://pypi.org/project/PaDEL-pywrapper/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/OlivierBeq/PaDEL_pywrapper/actions/workflows/tests.yml/badge.svg)](https://github.com/OlivierBeq/PaDEL_pywrapper/actions/workflows/tests.yml)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
<br>
</div>

A simple and reliable Python wrapper for calculating [PaDEL](https://doi.org/10.1002/jcc.21707) molecular descriptors and fingerprints. This library takes care of installing a matching Java runtime, dispatching molecules to the bundled PaDEL executable, and collecting the results into a tidy pandas DataFrame — so you can stay in RDKit/pandas-land.

## ✨ Features

- 🧬 **1875 descriptors & 12 fingerprint types** — 1444 1D/2D and 431 3D descriptors, plus FP, ExtFP, EStateFP, GraphFP, MACCSFP, PubchemFP, SubFP(C), KRFP(C) and AP2DFP(C) fingerprints, straight from PaDEL-Descriptor.
- ☕ **Zero Java setup** — automatically downloads, caches and reuses a matching JRE on first use; nothing to install by hand.
- ⚡ **Parallel by design** — spread the work across multiple CPU cores with configurable `njobs`/`chunksize`, each worker running its own single-core-pinned JVM.
- 🛡️ **Hang-proof** — some 3D descriptors (e.g. WHIM) can fail to converge on degenerate geometries; an optional per-descriptor `timeout` guarantees you get `NaN` back instead of a frozen process.
- 🧯 **Never silently misaligned** — molecules that fail, get skipped, or lack hydrogens/conformers are handled explicitly and reported, never dropped without a trace.
- 📊 **pandas-native output** — results come back as a ready-to-use `DataFrame`, one row per molecule.
- 🔍 **Rich metadata** — inspect each descriptor's description, sub-components and 3D requirement programmatically.
- 🧩 **Configurable fingerprints** — tune bit size and search depth for the fingerprints that support it.

## ✍️ Copyright and Citation Notice

Olivier J. M. Béquignon is **neither** the copyright holder of PaDEL **nor** responsible for it. The work carried out here concerns:

- the Python wrapper,
- the ePaDEL executable,
- the extendedlibpadeldescriptor library.

### Citing

If you use this wrapper in your research, please cite the original PaDEL publication in addition to this software package:

1. **Original PaDEL paper:**
   > Yap, C.W. (2011), PaDEL-descriptor: An open source software to calculate molecular descriptors and fingerprints. *Journal of Computational Chemistry*, 32(7), 1466–1474.
   > [DOI: 10.1002/jcc.21707](https://doi.org/10.1002/jcc.21707)

2. **This wrapper:**
   > Please refer to the repository at [github.com/OlivierBeq/PaDEL_pywrapper](https://github.com/OlivierBeq/PaDEL_pywrapper) for citation details.


## 📦 Installation


```bash
pip install padel-pywrapper
```

  Or from source:

```bash
git clone https://github.com/OlivierBeq/PaDEL_pywrapper.git
pip install ./PaDEL_pywrapper
```
## 🛠️ Requirements

- Python 3.11+
- [RDKit](https://www.rdkit.org/docs/Install.html)


## 💡 Usage

### 1D and 2D descriptors

  Descriptors of the module `PaDEL_pywrapper.descriptor` can be computed as follows:

```python
With pip:
from PaDEL_pywrapper import PaDEL
from PaDEL_pywrapper.descriptor import ALOGP, Crippen, FMF
from rdkit import Chem

smiles_list = [
# erlotinib
    "n1cnc(c2cc(c(cc12)OCCOC)OCCOC)Nc1cc(ccc1)C#C",
    # midecamycin
    "CCC(=O)O[C@@H]1CC(=O)O[C@@H](C/C=C/C=C/[C@@H]([C@@H](C[C@@H]([C@@H]([C@H]1OC)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C)O[C@H]3C[C@@]([C@H]([C@@H](O3)C)OC(=O)CC)(C)O)N(C)C)O)CC=O)C)O)C",
    # selenofolate
    "C1=CC(=CC=C1C(=O)NC(CCC(=O)OCC[Se]C#N)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N",
]
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

descriptors = [ALOGP, Crippen, FMF]

padel = PaDEL(descriptors)
print(padel.calculate(mols))
```

Instances of descriptors can be supplied as well:

```python
descriptors = [ALOGP(), Crippen(), FMF()]

padel = PaDEL(descriptors)
print(padel.calculate(mols))
```

To calculate all possible descriptors, import the `descriptors` list from the `PaDEL_pywrapper` module directly:

```python
from PaDEL_pywrapper import descriptors

padel = PaDEL(descriptors)
print(padel.calculate(mols))
```

### 3D descriptors

By default, the `ignore_3D` parameter is set to `True`, preventing any provided 3D descriptor from being calculated.

Should molecules with 3D coordinates be provided, one can turn on 3D descriptor calculation:

```python
from rdkit.Chem import AllChem
from PaDEL_pywrapper.descriptor import WHIM

mols = [Chem.AddHs(mol) for mol in mols]
_ = [AllChem.EmbedMolecule(mol) for mol in mols]

descriptors_3d = [ALOGP, Crippen, FMF, WHIM]

padel = PaDEL(descriptors_3d, ignore_3D=False)
print(padel.calculate(mols))
```

> ⚠️ A warning is raised if molecules lack hydrogens.
> ⚠️ An exception is raised if molecules do not have 3D coordinates — but only when a 3D
> descriptor (like `WHIM` above) is actually requested; a purely 2D descriptor list computes
> 2D coordinates on the fly instead.

```python
mol = Chem.MolFromSmiles('CCC')

padel = PaDEL(descriptors_3d, ignore_3D=False)
print(padel.calculate([mol]))
# ValueError: Cannot calculate descriptors for a conformer-less molecule
```

Some 3D descriptors (e.g. `WHIM`) rely on an eigenvalue decomposition that can fail to converge — and hang indefinitely — for molecules with degenerate 3D geometry (very few heavy atoms, near-linear/near-planar arrangements). 🛡️ Pass a `timeout` (in seconds) to guarantee that any descriptor exceeding it returns `NaN` for that molecule instead of blocking the whole batch:

```python
padel = PaDEL(descriptors_3d, ignore_3D=False)
print(padel.calculate(mols, timeout=20))
```

### Fingerprints

Fingerprints of the module `PaDEL_pywrapper.descriptor` can be computed as follows:

```python
from PaDEL_pywrapper.descriptor import GraphOnlyFP

fp = GraphOnlyFP

padel = PaDEL([fp], ignore_3D=False)
print(padel.calculate(mols))
```

Custom parameter sets can be provided for some fingerprints:

```python
fp = GraphOnlyFP(size=2048, searchDepth=8)

padel = PaDEL([fp], ignore_3D=False)
print(padel.calculate(mols))
```

### ⚡ Parallel processing

Speed things up by spreading molecules across several CPU cores. Each worker runs its own single-core-pinned JVM, so parallelism comes purely from the number of processes spawned — not from oversubscribing the host:

```python
padel = PaDEL(descriptors)
print(padel.calculate(mols, njobs=8))
```

By default, molecules are auto-balanced evenly across `njobs` workers (`chunksize=None`), which minimizes JVM startup overhead while keeping every worker busy — the fastest setting for most workloads. A fixed `chunksize` can be provided instead if finer control is needed.

### 🔍 Details about descriptors

Details about each descriptor and fingerprint can be obtained as follows:

```python
print(ALOGP.description)

print(GraphOnlyFP.description)
```

For full details about all descriptors, one can obtain the path to the original Excel file of the PaDEL descriptors with:

```python
print(padel.details)
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/OlivierBeq/PaDEL_pywrapper/blob/master/LICENSE) file for details.

## 📚 API Documentation

```python
def calculate(mols, show_banner=True, njobs=1, chunksize=None, timeout=None):
```

Calculates PaDEL molecular descriptors and/or fingerprints. Installs a matching JRE on first use if none is found.

#### Parameters

- ***mols  : Iterable[Chem.Mol]***
  RDKit molecule objects for which to obtain PaDEL descriptors.
- ***show_banner  : bool***
  Displays default notice about PaDEL descriptors.
- ***njobs  : int***
  Number of concurrent processes used to calculate descriptors in parallel; must not exceed the
  number of available CPU cores. Each spawned Java process is pinned to a single core
  (`-XX:ActiveProcessorCount=1`), since parallelism comes from spawning `njobs` OS processes
  rather than from letting each JVM oversubscribe the host's full core count.
- ***chunksize  : int | None***
  Number of molecules processed per worker process. If `None` (default), molecules are
  auto-balanced across all `njobs` workers so every worker gets work. Ignored if `njobs` is 1.
- ***timeout  : float | None***
  Maximum number of seconds allowed for each descriptor/fingerprint calculation for a single
  molecule. `None` (default) waits indefinitely. Guards against 3D descriptors (e.g. `WHIM`)
  that can hang on degenerate geometries — the offending molecule gets `NaN` instead of
  blocking the whole batch.
- ***return_type  : pd.DataFrame***
  Pandas DataFrame containing PaDEL molecular descriptors and/or fingerprints, one row per molecule.

________________

```python
Descriptor(name, is_3D)
```

Metadata holder for a single PaDEL descriptor family (e.g. `ALOGP`, `WHIM`). Instances are
pre-created and importable by name from `PaDEL_pywrapper.descriptor`.

#### Attributes

- ***name  : str***
  Name of the descriptor as known to PaDEL.
- ***is_3D  : bool***
  Whether the descriptor requires 3D molecular coordinates.
- ***subcomponents  : list[str]***
  Names of the individual output columns making up this descriptor.
- ***description  : pd.DataFrame***
  Human-readable description of each sub-component.

________________

```python
Fingerprint(name)
```

Metadata holder for a single PaDEL fingerprint type. Instances are pre-created and importable
by name from `PaDEL_pywrapper.descriptor`; configurable ones (`FP`, `GraphOnlyFP`) can be
called with `size`/`searchDepth` to override the defaults.

#### Attributes

- ***name  : str***
  Name of the fingerprint as known to PaDEL.
- ***size  : int | None***
  Number of bits, if configurable for this fingerprint type.
- ***searchDepth  : int | None***
  Search depth, if configurable for this fingerprint type.
