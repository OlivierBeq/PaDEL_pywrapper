[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# PaDEL Python wrapper

Python wrapper to ease the calculation of [PaDEL molecular descriptors](https://doi.org/10.1002/jcc.21707) and fingerprints.

## Copyright notice

Olivier J. M. Béquignon is **neither** the copyright holder of PaDEL **nor** responsible for it.

Only the Python wrapper is the work of Olivier J. M. Béquignon.

## Installation

From source:

    git clone https://github.com/OlivierBeq/PaDEL_pywrapper.git
    pip install ./PaDEL_pywrapper

with pip:

```bash
pip install padel-pywrapper
```

### Get started

#### 1D and 2D descriptors

Descriptors of the module `PaDEL_pywrapper.descriptors` can be computed as follows:

```python
from PaDEL_pywrapper import PaDEL
from PaDEL_pywrapper.descriptors import CDK_ALOGPDescriptor, CDK_CrippenDescriptor, CDK_FMFDescriptor
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

descriptors = [CDK_ALOGPDescriptor, CDK_CrippenDescriptor, CDK_FMFDescriptor]

padel = PaDEL(descriptors)
print(padel.calculate(mols))
```

Instances of descriptors can be supplied as well.

```python
descriptors = [CDK_ALOGPDescriptor(), CDK_CrippenDescriptor(), CDK_FMFDescriptor()]

padel = PaDEL(descriptors)
print(padel.calculate(mols))
```

To calculate all possible descriptors, import the `descriptors` list from the module `PaDEL_pywrapper` directly:


```python
from PaDEL_pywrapper import descriptors

padel = PaDEL(descriptors)
print(padel.calculate(mols))
```

#### 3D descriptors

By default, the `ignore_3D` parameter is set to `True`, preventing any provided 3D descriptor to be calculated.

Should molecules with 3D coordinates be provided, one can turn on 3D descriptor calculation.

```python
from rdkit.Chem import AllChem

mols = [Chem.AddHs(mol) for mol in mols]
_ = [AllChem.EmbedMolecule(mol) for mol in mols]

padel = PaDEL(descriptors, ignore_3D=False)
print(padel.calculate(mols))
```

:warning: An exception is raised if molecules do not have 3D coordinates.

```python
mol = Chem.MolFromSmiles('CCC')

padel = PaDEL(descriptors, ignore_3D=False)
print(padel.calculate([mol]))
# ValueError: Cannot calculate descriptors for a conformer-less molecule
```

#### Fingerprints


Fingerprints of the module `PaDEL_pywrapper.descriptors can be computed as follows:

```python
from PaDEL_pywrapper.descriptors import GraphOnlyFingerprint

fp = GraphOnlyFingerprint

padel = PaDEL([fp], ignore_3D=False)
print(padel.calculate(mols))
```

Custom parameter sets can be provided for some fingerprints:

```python
from PaDEL_pywrapper.descriptors import GraphOnlyFingerprint

fp = GraphOnlyFingerprint()
fp.set_params({'size': 2048, 'searchDepth': 8})

padel = PaDEL([fp], ignore_3D=False)
print(padel.calculate(mols))
```

### Other parameters

```python
class PaDEL:
    ...
    def calculate(self, mols: Iterable[Chem.Mol], show_banner: bool = True, njobs: int = 1, chunksize: int = 100):
```

#### Parameters

- ***mols  : Iterable[Chem.Mol]***  
  RDKit molecule objects for which to obtain PaDEL descriptors.
- ***show_banner  : bool***  
  Displays default notice about PaDEL descriptors.
- ***njobs  : int***  
  Maximum number of simultaneous processes. Ignored if `self.descriptors` are instances and not class names.
- ***chunksize  : int***  
  Maximum number of molecules each process is charged of. Ignored if `self.descriptors` are instances and not class names.

### Details about descriptors

The path to the spreadsheet containing all details about the PaDEL descriptors can be found using:

```python
print(padel.details)
```
