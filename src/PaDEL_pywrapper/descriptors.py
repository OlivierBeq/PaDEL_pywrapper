# -*- coding: utf-8

"""Available PaDEL descriptors."""

import os
import re
import sys
from abc import ABC, abstractmethod
from typing import Dict, List, Union

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from jpype import isJVMStarted, startJVM, getDefaultJVMPath, JPackage, types as jtypes

from .utils import parse_numeric, parse_numeric_to_null


class Descriptor(ABC):
    """Primitive to allow multiprocessed calculation of PaDEL descriptors."""

    def __init__(self, is3D: bool = False):
        """Instantiate a new Descriptor

        :param is3D: is a 3D descriptor
        """
        self.is3D = is3D
        self._load_resources()

    @property
    @abstractmethod
    def descriptor(self):
        pass

    def _load_resources(self):
        """Load all resources required for the descriptor calculation."""
        if not isJVMStarted():
            # Path to the PaDEL java library
            padel_path = os.path.join(__file__, os.pardir, 'PaDEL-Descriptor', 'lib', 'libPaDEL-Descriptor.jar')
            # Path to the additions to the PaDEL java library
            epadel_path = os.path.join(__file__, os.pardir, 'PaDEL-Descriptor', 'lib', 'elibPaDEL-Descriptor.jar')
            # Platform dependent java classpath separator
            separator = ':' if sys.platform == 'linux' else ';'
            # Start JVM
            startJVM(getDefaultJVMPath(), "-ea", f"-Djava.class.path={padel_path}{separator}{epadel_path}")
        # Java modules to be accessed
        self.padel = JPackage('libpadeldescriptor')
        self.epadel = JPackage('extendedlibpadeldescriptor')
        self.cdk = JPackage('org.openscience.cdk')
        self.io = JPackage('java.io')
        self.lang = JPackage('java.lang')

    def _prepare_mols(self, mols: List[Chem.Mol]):
        """Convert RDKit molecules into CDK molecules."""
        cdk_mols = []
        for mol in mols:
            confs = list(mol.GetConformers())
            # If mols are not 3D, compute 2D coords
            if not (len(confs) > 0 and confs[-1].Is3D()):
                # If descriptor is 3D raise error
                if self.is3D:
                    raise ValueError('Cannot calculate descriptors for a conformer-less molecule')
                AllChem.Compute2DCoords(mol)

            # Pass to CDK through SDF parser
            buffer = self.io.ByteArrayInputStream(self.lang.String(Chem.MolToV3KMolBlock(mol)).getBytes("UTF-8"))
            parser = self.cdk.io.MDLV3000Reader(buffer)
            cdk_mol = parser.read(self.cdk.Molecule())
            cdk_mols.append(cdk_mol)
        return cdk_mols

    def calculate(self, mol: Union[Chem.Mol, List[Chem.Mol]]):
        """Obtain values of the descriptors.

        :param mol: molecule or molecules to obtain the descriptor of
        """
        # Load resources
        self._load_resources()
        # Ensure input is a list
        if isinstance(mol, Chem.Mol):
            mol = [mol]
        # Convert to CDK
        cdk_mols = self._prepare_mols(mol)

        # Ensure empty fallback is ready
        self._empty_val()

        # Calculate values
        data = []
        descriptor_names = []
        for i, cdk_mol in enumerate(cdk_mols):
            mol_descs = []
            # CDK descriptor
            if hasattr(self.descriptor, 'setMolecule'):
                self.descriptor.setMolecule(self.cdk.Molecule())
                self.descriptor.setMolecule(cdk_mol)
                try:
                    self.descriptor.run()
                    # Get descriptor names for 1st molecule only
                    if i == 0:
                        descriptor_names.extend(str(x) for x in self.descriptor.getDescriptorNames())
                    values = list(parse_numeric(str(x)) for x in self.descriptor.getDescriptorValues())
                    mol_descs.extend(values)
                except:
                    if i == 0:
                        descriptor_names.extend(self._names)
                    mol_descs.extend(self._empty)
            # PaDEL custom descriptor
            elif hasattr(self.descriptor, 'calculate'):
                res = self.descriptor.calculate(cdk_mol)
                res_value = res.getValue()
                # Get descriptor names for 1st molecule only
                if i == 0:
                    descriptor_names.extend(str(x) for x in res.getNames())
                values = list(parse_numeric(str(res_value.get(i))) for i in range(res_value.length()))
                mol_descs.extend(values)
            else:
                print('FAILED: ', type(self.descriptor).__name__)
            data.append(mol_descs)

        # Format as pandas DataFrame
        return pd.DataFrame(data, columns=descriptor_names)

    @property
    def name(self) -> str:
        """Obtain the name of a descriptor."""
        return re.sub('^e?CDK_|Descripotr$', '',
                      type(self).__name__.split('.')[-1]).replace('Fingerprinter', 'Fingerprint')

    def _empty_val(self):
        """Return an empty row for the current descriptor."""
        if not hasattr(self, '_empty'):
            # Use proxy molecule
            if self.is3D:
                mol = Chem.MolFromSmiles('CC')
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol)
                cdk_mol = self._prepare_mols([mol])[0]
            else:
                cdk_mol = self._prepare_mols([Chem.MolFromSmiles('C')])[0]
            if hasattr(self.descriptor, 'setMolecule'):
                self.descriptor.setMolecule(self.cdk.Molecule())
                self.descriptor.setMolecule(cdk_mol)
                self.descriptor.run()
                # Get descriptor names for 1st molecule only
                self._names = list(map(str, self.descriptor.getDescriptorNames()))
                self._empty = list(map(parse_numeric_to_null, self.descriptor.getDescriptorValues()))
            # PaDEL custom descriptor
            elif hasattr(self.descriptor, 'calculate'):
                res = self.descriptor.calculate(cdk_mol)
                res_value = res.getValue()
                # Get descriptor names for 1st molecule only
                self._names = list(map(str, res.getNames()))
                self._empty = list(parse_numeric_to_null(str(res_value.get(i))) for i in range(res_value.length()))


class CDK_AcidicGroupCountDescriptor(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AcidicGroupCountDescriptor()
        return self._descriptor


class CDK_ALOGPDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_ALOGPDescriptor()
        return self._descriptor


class CDK_APolDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_APolDescriptor()
        return self._descriptor


class CDK_AromaticAtomsCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AromaticAtomsCountDescriptor()
        return self._descriptor


class CDK_AromaticBondsCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AromaticBondsCountDescriptor()
        return self._descriptor


class CDK_AllAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AtomCountDescriptor(jtypes.JArray(str, 1)(['*']))
        return self._descriptor


class CDK_HeavyAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_HeavyAtomCountDescriptor()
        return self._descriptor


class CDK_HAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AtomCountDescriptor(jtypes.JArray(str, 1)(['H']))
        return self._descriptor


class CDK_BAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AtomCountDescriptor(jtypes.JArray(str, 1)(['B']))
        return self._descriptor


class CDK_CAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AtomCountDescriptor(jtypes.JArray(str, 1)(['C']))
        return self._descriptor


class CDK_NAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AtomCountDescriptor(jtypes.JArray(str, 1)(['N']))
        return self._descriptor


class CDK_SAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AtomCountDescriptor(jtypes.JArray(str, 1)(['S']))
        return self._descriptor


class CDK_OAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AtomCountDescriptor(jtypes.JArray(str, 1)(['O']))
        return self._descriptor


class CDK_PAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AtomCountDescriptor(jtypes.JArray(str, 1)(['P']))
        return self._descriptor


class CDK_FAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AtomCountDescriptor(jtypes.JArray(str, 1)(['F']))
        return self._descriptor


class CDK_ClAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AtomCountDescriptor(jtypes.JArray(str, 1)(['Cl']))
        return self._descriptor


class CDK_BrAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AtomCountDescriptor(jtypes.JArray(str, 1)(['Br']))
        return self._descriptor


class CDK_IAtomCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AtomCountDescriptor(jtypes.JArray(str, 1)(['I']))
        return self._descriptor


class CDK_HalogenCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_HalogenCountDescriptor()
        return self._descriptor


class CDK_AutocorrelationDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_AutocorrelationDescriptor()
        return self._descriptor


class BaryszMatrixDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.BaryszMatrixDescriptor()
        return self._descriptor


class CDK_BasicGroupCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_BasicGroupCountDescriptor()
        return self._descriptor


class CDK_BCUTDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_BCUTDescriptor()
        return self._descriptor


class CDK_BondCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_BondCountDescriptor()
        return self._descriptor


class CDK_BPolDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_BPolDescriptor()
        return self._descriptor


class CDK_BurdenModifiedEigenvaluesDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_BurdenModifiedEigenvaluesDescriptor()
        return self._descriptor


class CDK_CarbonTypesDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_CarbonTypesDescriptor()
        return self._descriptor


class CDK_ChiChainDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_ChiChainDescriptor()
        return self._descriptor


class CDK_ChiClusterDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_ChiClusterDescriptor()
        return self._descriptor


class CDK_ChiPathClusterDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_ChiPathClusterDescriptor()
        return self._descriptor


class CDK_ChiPathDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_ChiPathDescriptor()
        return self._descriptor


class CDK_ConstitutionalDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_ConstitutionalDescriptor()
        return self._descriptor


class CDK_CrippenDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_CrippenDescriptor()
        return self._descriptor


class CDK_DetourMatrixDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_DetourMatrixDescriptor()
        return self._descriptor


class CDK_EccentricConnectivityIndexDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_EccentricConnectivityIndexDescriptor()
        return self._descriptor


class CDK_EStateAtomTypeDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_EStateAtomTypeDescriptor()
        return self._descriptor


class CDK_ExtendedTopochemicalAtomDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_ExtendedTopochemicalAtomDescriptor()
        return self._descriptor


class CDK_FMFDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_FMFDescriptor()
        return self._descriptor


class CDK_FragmentComplexityDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_FragmentComplexityDescriptor()
        return self._descriptor


class CDK_HBondAcceptorCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_HBondAcceptorCountDescriptor()
        return self._descriptor


class CDK_HBondDonorCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_HBondDonorCountDescriptor()
        return self._descriptor


class CDK_HybridizationRatioDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_HybridizationRatioDescriptor()
        return self._descriptor


class CDK_InformationContentDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_InformationContentDescriptor()
        return self._descriptor


class CDK_KappaShapeIndicesDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_KappaShapeIndicesDescriptor()
        return self._descriptor


class CDK_LargestChainDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_LargestChainDescriptor()
        return self._descriptor


class CDK_LargestPiSystemDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_LargestPiSystemDescriptor()
        return self._descriptor


class CDK_LongestAliphaticChainDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_LongestAliphaticChainDescriptor()
        return self._descriptor


class CDK_MannholdLogPDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_MannholdLogPDescriptor()
        return self._descriptor


class CDK_McGowanVolumeDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_McGowanVolumeDescriptor()
        return self._descriptor


class CDK_MDEDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_MDEDescriptor()
        return self._descriptor


class CDK_MLFERDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_MLFERDescriptor()
        return self._descriptor


class CDK_PathCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_PathCountDescriptor()
        return self._descriptor


class CDK_PetitjeanNumberDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_PetitjeanNumberDescriptor()
        return self._descriptor


class CDK_RingCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_RingCountDescriptor()
        return self._descriptor


class CDK_RotatableBondsCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_RotatableBondsCountDescriptor()
        return self._descriptor


class CDK_RuleOfFiveDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_RuleOfFiveDescriptor()
        return self._descriptor


class CDK_TopologicalDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_TopologicalDescriptor()
        return self._descriptor


class CDK_TopologicalChargeDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_TopologicalChargeDescriptor()
        return self._descriptor


class CDK_TopologicalDistanceMatrixDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_TopologicalDistanceMatrixDescriptor()
        return self._descriptor


class CDK_TPSADescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_TPSADescriptor()
        return self._descriptor


class CDK_VABCDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_VABCDescriptor()
        return self._descriptor


class CDK_VAdjMaDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_VAdjMaDescriptor()
        return self._descriptor


class CDK_WalkCountDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_WalkCountDescriptor()
        return self._descriptor


class CDK_WeightDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_WeightDescriptor()
        return self._descriptor


class CDK_WeightedPathDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_WeightedPathDescriptor()
        return self._descriptor


class CDK_WienerNumbersDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_WienerNumbersDescriptor()
        return self._descriptor


class CDK_XLogPDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_XLogPDescriptor(jtypes.JArray(jtypes.JObject, 1)([False, True]))
        return self._descriptor


class CDK_ZagrebIndexDescriptor(Descriptor):
    def __init__(self):
        super().__init__(False)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_ZagrebIndexDescriptor()
        return self._descriptor


class CDK_Autocorrelation3DDescriptor(Descriptor):
    def __init__(self):
        super().__init__(True)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_Autocorrelation3DDescriptor()
        return self._descriptor


class CDK_CPSADescriptor(Descriptor):
    def __init__(self):
        super().__init__(True)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_CPSADescriptor()
        return self._descriptor


class CDK_GravitationalIndexDescriptor(Descriptor):
    def __init__(self):
        super().__init__(True)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_GravitationalIndexDescriptor()
        return self._descriptor


class CDK_LengthOverBreadthDescriptor(Descriptor):
    def __init__(self):
        super().__init__(True)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_LengthOverBreadthDescriptor()
        return self._descriptor


class CDK_MomentOfInertiaDescriptor(Descriptor):
    def __init__(self):
        super().__init__(True)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_MomentOfInertiaDescriptor()
        return self._descriptor


class CDK_PetitjeanShapeIndexDescriptor(Descriptor):
    def __init__(self):
        super().__init__(True)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_PetitjeanShapeIndexDescriptor()
        return self._descriptor


class CDK_RDFDescriptor(Descriptor):
    def __init__(self):
        super().__init__(True)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_RDFDescriptor()
        return self._descriptor


class CDK_WHIMDescriptor(Descriptor):
    def __init__(self):
        super().__init__(True)

    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.padel.CDK_WHIMDescriptor()
        return self._descriptor


class eCDK_Fingerprinter(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.epadel.eCDK_Fingerprinter()
        return self._descriptor

    def set_params(self, params: Dict[str, int] = {'size': 1024, 'searchDepth': 7}):
        default = {'size': 1024, 'searchDepth': 7}
        default.update(params)
        self._descriptor = self.epadel.eCDK_Fingerprinter(default['size'], default['searchDepth'])


class eCDK_ExtendedFingerprinter(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.epadel.eCDK_ExtendedFingerprinter()
        return self._descriptor

    def set_params(self, params: Dict[str, int] = {'size': 1024, 'searchDepth': 7}):
        default = {'size': 1024, 'searchDepth': 7}
        default.update(params)
        self._descriptor = self.epadel.eCDK_ExtendedFingerprinter(default['size'], default['searchDepth'])


class eCDK_EStateFingerprinter(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.epadel.eCDK_EStateFingerprinter()
        return self._descriptor


class eCDK_GraphOnlyFingerprinter(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.epadel.eCDK_GraphOnlyFingerprinter()
        return self._descriptor

    def set_params(self, params: Dict[str, int] = {'size': 1024, 'searchDepth': 7}):
        default = {'size': 1024, 'searchDepth': 7}
        default.update(params)
        self._descriptor = self.epadel.eCDK_GraphOnlyFingerprinter(default['size'], default['searchDepth'])


class eCDK_MACCSFingerprinter(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.epadel.eCDK_MACCSFingerprinter()
        return self._descriptor


class eCDK_PubchemFingerprinter(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.epadel.eCDK_PubchemFingerprinter()
        return self._descriptor


class eCDK_SubstructureFingerprinter(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.epadel.eCDK_SubstructureFingerprinter()
        return self._descriptor


class eCDK_KlekotaRothFingerprinter(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.epadel.eCDK_KlekotaRothFingerprinter()
        return self._descriptor


class eCDK_AtomPairs2DFingerprinter(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.epadel.eCDK_AtomPairs2DFingerprinter()
        return self._descriptor


class eCDK_SubstructureFingerprintCount(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.epadel.eCDK_SubstructureFingerprintCount()
        return self._descriptor


class eCDK_KlekotaRothFingerprintCount(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.epadel.eCDK_KlekotaRothFingerprintCount()
        return self._descriptor


class eCDK_AtomPairs2DFingerprintCount(Descriptor):
    @property
    def descriptor(self):
        if not hasattr(self, '_descriptor'):
            self._descriptor = self.epadel.eCDK_AtomPairs2DFingerprintCount()
        return self._descriptor


_descs_2D = [CDK_AcidicGroupCountDescriptor,
             CDK_ALOGPDescriptor,
             CDK_APolDescriptor,
             CDK_AromaticAtomsCountDescriptor,
             CDK_AromaticBondsCountDescriptor,
             CDK_AllAtomCountDescriptor,
             CDK_HeavyAtomCountDescriptor,
             CDK_HAtomCountDescriptor,
             CDK_BAtomCountDescriptor,
             CDK_CAtomCountDescriptor,
             CDK_NAtomCountDescriptor,
             CDK_OAtomCountDescriptor,
             CDK_SAtomCountDescriptor,
             CDK_PAtomCountDescriptor,
             CDK_HalogenCountDescriptor,
             CDK_AutocorrelationDescriptor,
             BaryszMatrixDescriptor,
             CDK_BasicGroupCountDescriptor,
             CDK_BCUTDescriptor,
             CDK_BondCountDescriptor,
             CDK_BPolDescriptor,
             CDK_BurdenModifiedEigenvaluesDescriptor,
             CDK_CarbonTypesDescriptor,
             CDK_ChiChainDescriptor,
             CDK_ChiClusterDescriptor,
             CDK_ChiPathClusterDescriptor,
             CDK_ChiPathDescriptor,
             CDK_ConstitutionalDescriptor,
             CDK_CrippenDescriptor,
             CDK_DetourMatrixDescriptor,
             CDK_EccentricConnectivityIndexDescriptor,
             CDK_EStateAtomTypeDescriptor,
             CDK_ExtendedTopochemicalAtomDescriptor,
             CDK_FMFDescriptor,
             CDK_FragmentComplexityDescriptor,
             CDK_HBondAcceptorCountDescriptor,
             CDK_HBondDonorCountDescriptor,
             CDK_HybridizationRatioDescriptor,
             CDK_InformationContentDescriptor,
             CDK_KappaShapeIndicesDescriptor,
             CDK_LargestChainDescriptor,
             CDK_LargestPiSystemDescriptor,
             CDK_LongestAliphaticChainDescriptor,
             CDK_MannholdLogPDescriptor,
             CDK_McGowanVolumeDescriptor,
             CDK_MDEDescriptor,
             CDK_MLFERDescriptor,
             CDK_PathCountDescriptor,
             CDK_PetitjeanNumberDescriptor,
             CDK_RingCountDescriptor,
             CDK_RotatableBondsCountDescriptor,
             CDK_RuleOfFiveDescriptor,
             CDK_TopologicalDescriptor,
             CDK_TopologicalChargeDescriptor,
             CDK_TopologicalDistanceMatrixDescriptor,
             CDK_TPSADescriptor,
             CDK_VABCDescriptor,
             CDK_VAdjMaDescriptor,
             CDK_WalkCountDescriptor,
             CDK_WeightDescriptor,
             CDK_WeightedPathDescriptor,
             CDK_WienerNumbersDescriptor,
             CDK_XLogPDescriptor,
             CDK_ZagrebIndexDescriptor]

_descs_3D = [CDK_Autocorrelation3DDescriptor,
             CDK_CPSADescriptor,
             CDK_GravitationalIndexDescriptor,
             CDK_LengthOverBreadthDescriptor,
             CDK_MomentOfInertiaDescriptor,
             CDK_PetitjeanShapeIndexDescriptor,
             CDK_RDFDescriptor,
             CDK_WHIMDescriptor
             ]

descriptors = _descs_2D + _descs_3D

_fingerprints = [eCDK_Fingerprinter,
                 eCDK_ExtendedFingerprinter,
                 eCDK_EStateFingerprinter,
                 eCDK_GraphOnlyFingerprinter,
                 eCDK_MACCSFingerprinter,
                 eCDK_PubchemFingerprinter,
                 eCDK_SubstructureFingerprinter,
                 eCDK_KlekotaRothFingerprinter,
                 eCDK_AtomPairs2DFingerprinter,
                 eCDK_SubstructureFingerprintCount,
                 eCDK_KlekotaRothFingerprintCount,
                 eCDK_AtomPairs2DFingerprintCount,
                 ]

Fingerprint = _fingerprints[0]
ExtendedFingerprint = _fingerprints[1]
EStateFingerprint = _fingerprints[2]
GraphOnlyFingerprint = _fingerprints[3]
MACCSFingerprint = _fingerprints[4]
PubchemFingerprint = _fingerprints[5]
SubstructureFingerprint = _fingerprints[6]
KlekotaRothFingerprint = _fingerprints[7]
AtomPairs2DFingerprint = _fingerprints[8]
SubstructureFingerprintCount = _fingerprints[9]
KlekotaRothFingerprintCount = _fingerprints[10]
AtomPairs2DFingerprintCount = _fingerprints[11]
