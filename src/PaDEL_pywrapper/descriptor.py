
# -*- coding: utf-8

"""Python wrapper for PaDEL descriptors"""

import os

import numpy as np
import pandas as pd


class Descriptor:
    """Primitive PaDEL descriptor."""

    def __init__(self, name: str, is_3D: bool):
        """Instantiate a descriptor.

        :param name: Name of the descriptor
        :param is_3D: Is the descriptor requiring 3D molecular coordinates
        :param type: Type of descriptor
        :param description: Description of the descriptor
        """
        self.name = name
        self.is_3D = is_3D
        # Get names of bits (e.g. for AtomPairs2DFP: AP2DFP)
        self._get_sub_descriptor_names('padel_descs.tsv')

    def _get_sub_descriptor_names(self, path: str):
        """Load names of the subcomponents making this descriptor."""
        path = os.path.abspath(os.path.join(os.path.dirname(__file__), path))
        descs = pd.read_csv(path, sep='\t')
        self.subcomponents = descs[descs.descriptor == self.name].name.tolist()
        self.description = descs[descs.descriptor == self.name][['name', 'description']]


class Fingerprint:
    """Primitive PaDEL fingerprint."""

    def __init__(self, name: str):
        """Instantiate a fingerprint.

        :param name: Name of the fingerprint
        :param type: Name ePaDEL uses to compute the fingerprint
        """
        self.name = name
        self.is_3D = False
        # Get names of bits (e.g. for Atom: CrippenLogP & CrippenMR)
        self._get_sub_descriptor_names('padel_fps.tsv')

    def _get_sub_descriptor_names(self, path: str):
        """Load names of the subcomponents making this descriptor."""
        path = os.path.abspath(os.path.join(os.path.dirname(__file__), path))
        fps = pd.read_csv(path, sep='\t')
        self.bit_prefix = fps[fps.fingerprint == self.name].name.values[0]
        nbits = fps[fps.fingerprint == self.name]['fixed size'].values[0]
        self.nBits = nbits if nbits != '-' else np.NaN
        self.short_name = fps[fps.fingerprint == self.name]['name'].item()
        self.description = fps[fps.fingerprint == self.name]['description']

    def __call__(self, size: int = 1024, searchDepth: int = 7):
        """Define size and search depth for CDK and Graph only fingerprints.

        :param size: Number of bits in the CDK or Graph only fingerprints (ignored for others)
        :param searchDepth: Search depth for the CDK or Graph only fingerprints (ignored for others)
        """
        fp = Fingerprint(self.name)
        if self.name in ['CDK fingerprint', 'CDK graph only fingerprint']:
            fp.size = fp.nBits = size
            fp.searchDepth = searchDepth
        return fp


AcidicGroupCount = Descriptor('AcidicGroupCount', False)
ALOGP = Descriptor('ALOGP', False)
APol = Descriptor('APol', False)
AromaticAtomsCount = Descriptor('AromaticAtomsCount', False)
AromaticBondsCount = Descriptor('AromaticBondsCount', False)
AtomCount = Descriptor('AtomCount', False)
Autocorrelation = Descriptor('Autocorrelation', False)
BaryszMatrix = Descriptor('BaryszMatrix', False)
BasicGroupCount = Descriptor('BasicGroupCount', False)
BCUT = Descriptor('BCUT', False)
BondCount = Descriptor('BondCount', False)
BPol = Descriptor('BPol', False)
BurdenModifiedEigenvalues = Descriptor('BurdenModifiedEigenvalues', False)
CarbonTypes = Descriptor('CarbonTypes', False)
ChiChain = Descriptor('ChiChain', False)
ChiCluster = Descriptor('ChiCluster', False)
ChiPathCluster = Descriptor('ChiPathCluster', False)
ChiPath = Descriptor('ChiPath', False)
Constitutional = Descriptor('Constitutional', False)
Crippen = Descriptor('Crippen', False)
DetourMatrix = Descriptor('DetourMatrix', False)
EccentricConnectivityIndex = Descriptor('EccentricConnectivityIndex', False)
EStateAtomType = Descriptor('ElectrotopologicalStateAtomType', False)
ExtendedTopochemicalAtom = Descriptor('ExtendedTopochemicalAtom', False)
FMF = Descriptor('FMF', False)
FragmentComplexity = Descriptor('FragmentComplexity', False)
HBondAcceptorCount = Descriptor('HBondAcceptorCount', False)
HBondDonorCount = Descriptor('HBondDonorCount', False)
HybridizationRatio = Descriptor('HybridizationRatio', False)
InformationContent = Descriptor('InformationContent', False)
KappaShapeIndices = Descriptor('KappaShapeIndices', False)
LargestChain = Descriptor('LargestChain', False)
LargestPiSystem = Descriptor('LargestPiSystem', False)
LongestAliphaticChain = Descriptor('LongestAliphaticChain', False)
MannholdLogP = Descriptor('MannholdLogP', False)
McGowanVolume = Descriptor('McGowanVolume', False)
MDE = Descriptor('MDE', False)
MLFER = Descriptor('MLFER', False)
PathCount = Descriptor('PathCount', False)
PetitjeanNumber = Descriptor('PetitjeanNumber', False)
RingCount = Descriptor('RingCount', False)
RotatableBondsCount = Descriptor('RotatableBondsCount', False)
RuleOfFive = Descriptor('RuleOfFive', False)
Topological = Descriptor('Topological', False)
TopologicalCharge = Descriptor('TopologicalCharge', False)
TopologicalDistanceMatrix = Descriptor('TopologicalDistanceMatrix', False)
TPSA = Descriptor('TPSA', False)
VABC = Descriptor('VABC', False)
VAdjMa = Descriptor('VAdjMa', False)
WalkCount = Descriptor('WalkCount', False)
Weight = Descriptor('Weight', False)
WeightedPath = Descriptor('WeightedPath', False)
WienerNumbers = Descriptor('WienerNumbers', False)
XLogP = Descriptor('XLogP', False)
ZagrebIndex = Descriptor('ZagrebIndex', False)
Autocorrelation3D = Descriptor('Autocorrelation3D', True)
CPSA = Descriptor('CPSA', True)
GravitationalIndex = Descriptor('GravitationalIndex', True)
LengthOverBreadth = Descriptor('LengthOverBreadth', True)
MomentOfInertia = Descriptor('MomentOfInertia', True)
PetitjeanShapeIndex = Descriptor('PetitjeanShapeIndex', True)
RDF = Descriptor('RDF', True)
WHIM = Descriptor('WHIM', True)
FP = Fingerprint('CDK fingerprint')
ExtendedFP = Fingerprint('CDK extended fingerprint')
EStateFP = Fingerprint('Estate fingerprint')
GraphOnlyFP = Fingerprint('CDK graph only fingerprint')
MACCSFP = Fingerprint('MACCS fingerprint')
PubchemFP = Fingerprint('Pubchem fingerprint')
SubstructureFP = Fingerprint('Substructure fingerprint')
SubstructureFPCount = Fingerprint('Substructure fingerprint count')
KlekotaRothFP = Fingerprint('Klekota-Roth fingerprint')
KlekotaRothFPCount = Fingerprint('Klekota-Roth fingerprint count')
AtomPairs2DFP = Fingerprint('2D atom pairs fingerprint')
AtomPairs2DFPCount = Fingerprint('2D atom pairs fingerprint count')


descriptors = [AcidicGroupCount, ALOGP, APol, AromaticAtomsCount, AromaticBondsCount, AtomCount,
               Autocorrelation, BaryszMatrix, BasicGroupCount, BCUT, BondCount, BPol,
               BurdenModifiedEigenvalues, CarbonTypes, ChiChain, ChiCluster, ChiPathCluster, ChiPath,
               Constitutional, Crippen, DetourMatrix, EccentricConnectivityIndex, EStateAtomType,
               ExtendedTopochemicalAtom, FMF, FragmentComplexity, HBondAcceptorCount, HBondDonorCount,
               HybridizationRatio, InformationContent, KappaShapeIndices, LargestChain, LargestPiSystem,
               LongestAliphaticChain, MannholdLogP, McGowanVolume, MDE, MLFER, PathCount, PetitjeanNumber,
               RingCount, RotatableBondsCount, RuleOfFive, Topological, TopologicalCharge,
               TopologicalDistanceMatrix, TPSA, VABC, VAdjMa, WalkCount, Weight, WeightedPath, WienerNumbers,
               XLogP, ZagrebIndex, Autocorrelation3D, CPSA, GravitationalIndex, LengthOverBreadth,
               MomentOfInertia, PetitjeanShapeIndex, RDF, WHIM]

_fingerprints = [FP, ExtendedFP, EStateFP, GraphOnlyFP, MACCSFP, PubchemFP, SubstructureFP,
                 SubstructureFPCount, KlekotaRothFP, KlekotaRothFPCount, AtomPairs2DFP, AtomPairs2DFPCount]
