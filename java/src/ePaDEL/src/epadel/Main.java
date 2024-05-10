
package epadel;

import java.io.FileInputStream;
import java.io.IOException;
import java.lang.String;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.lang.ArrayIndexOutOfBoundsException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.smiles.SmilesParser;

import libpadeldescriptor.CDK_Descriptor;
import libpadeldescriptor.CDK_AcidicGroupCountDescriptor;
import libpadeldescriptor.CDK_ALOGPDescriptor;
import libpadeldescriptor.CDK_APolDescriptor;
import libpadeldescriptor.CDK_AromaticAtomsCountDescriptor;
import libpadeldescriptor.CDK_AromaticBondsCountDescriptor;
import libpadeldescriptor.CDK_AtomCountDescriptor;
import libpadeldescriptor.CDK_HeavyAtomCountDescriptor;
import libpadeldescriptor.CDK_HalogenCountDescriptor;
import libpadeldescriptor.CDK_AutocorrelationDescriptor;
import libpadeldescriptor.CDK_BaryszMatrixDescriptor;
import libpadeldescriptor.CDK_BasicGroupCountDescriptor;
import libpadeldescriptor.CDK_BCUTDescriptor;
import libpadeldescriptor.CDK_BondCountDescriptor;
import libpadeldescriptor.CDK_BPolDescriptor;
import libpadeldescriptor.CDK_BurdenModifiedEigenvaluesDescriptor;
import libpadeldescriptor.CDK_CarbonTypesDescriptor;
import libpadeldescriptor.CDK_ChiChainDescriptor;
import libpadeldescriptor.CDK_ChiClusterDescriptor;
import libpadeldescriptor.CDK_ChiPathClusterDescriptor;
import libpadeldescriptor.CDK_ChiPathDescriptor;
import libpadeldescriptor.CDK_ConstitutionalDescriptor;
import libpadeldescriptor.CDK_CrippenDescriptor;
import libpadeldescriptor.CDK_DetourMatrixDescriptor;
import libpadeldescriptor.CDK_EccentricConnectivityIndexDescriptor;
import libpadeldescriptor.CDK_EStateAtomTypeDescriptor;
import libpadeldescriptor.CDK_ExtendedTopochemicalAtomDescriptor;
import libpadeldescriptor.CDK_FMFDescriptor;
import libpadeldescriptor.CDK_FragmentComplexityDescriptor;
import libpadeldescriptor.CDK_HBondAcceptorCountDescriptor;
import libpadeldescriptor.CDK_HBondDonorCountDescriptor;
import libpadeldescriptor.CDK_HybridizationRatioDescriptor;
import libpadeldescriptor.CDK_InformationContentDescriptor;
import libpadeldescriptor.CDK_KappaShapeIndicesDescriptor;
import libpadeldescriptor.CDK_LargestChainDescriptor;
import libpadeldescriptor.CDK_LargestPiSystemDescriptor;
import libpadeldescriptor.CDK_LongestAliphaticChainDescriptor;
import libpadeldescriptor.CDK_MannholdLogPDescriptor;
import libpadeldescriptor.CDK_McGowanVolumeDescriptor;
import libpadeldescriptor.CDK_MDEDescriptor;
import libpadeldescriptor.CDK_MLFERDescriptor;
import libpadeldescriptor.CDK_PathCountDescriptor;
import libpadeldescriptor.CDK_PetitjeanNumberDescriptor;
import libpadeldescriptor.CDK_RingCountDescriptor;
import libpadeldescriptor.CDK_RotatableBondsCountDescriptor;
import libpadeldescriptor.CDK_RuleOfFiveDescriptor;
import libpadeldescriptor.CDK_TopologicalDescriptor;
import libpadeldescriptor.CDK_TopologicalChargeDescriptor;
import libpadeldescriptor.CDK_TopologicalDistanceMatrixDescriptor;
import libpadeldescriptor.CDK_TPSADescriptor;
import libpadeldescriptor.CDK_VABCDescriptor;
import libpadeldescriptor.CDK_VAdjMaDescriptor;
import libpadeldescriptor.CDK_WalkCountDescriptor;
import libpadeldescriptor.CDK_WeightDescriptor;
import libpadeldescriptor.CDK_WeightedPathDescriptor;
import libpadeldescriptor.CDK_WienerNumbersDescriptor;
import libpadeldescriptor.CDK_XLogPDescriptor;
import libpadeldescriptor.CDK_ZagrebIndexDescriptor;
import libpadeldescriptor.CDK_Autocorrelation3DDescriptor;
import libpadeldescriptor.CDK_CPSADescriptor;
import libpadeldescriptor.CDK_GravitationalIndexDescriptor;
import libpadeldescriptor.CDK_LengthOverBreadthDescriptor;
import libpadeldescriptor.CDK_MomentOfInertiaDescriptor;
import libpadeldescriptor.CDK_PetitjeanShapeIndexDescriptor;
import libpadeldescriptor.CDK_RDFDescriptor;
import libpadeldescriptor.CDK_WHIMDescriptor;
import extendedlibpadeldescriptor.eCDK_IFingerprint;
import extendedlibpadeldescriptor.eCDK_Fingerprinter;
import extendedlibpadeldescriptor.eCDK_ExtendedFingerprinter;
import extendedlibpadeldescriptor.eCDK_EStateFingerprinter;
import extendedlibpadeldescriptor.eCDK_GraphOnlyFingerprinter;
import extendedlibpadeldescriptor.eCDK_MACCSFingerprinter;
import extendedlibpadeldescriptor.eCDK_PubchemFingerprinter;
import extendedlibpadeldescriptor.eCDK_SubstructureFingerprinter;
import extendedlibpadeldescriptor.eCDK_KlekotaRothFingerprinter;
import extendedlibpadeldescriptor.eCDK_AtomPairs2DFingerprinter;
import extendedlibpadeldescriptor.eCDK_SubstructureFingerprintCount;
import extendedlibpadeldescriptor.eCDK_KlekotaRothFingerprintCount;
import extendedlibpadeldescriptor.eCDK_AtomPairs2DFingerprintCount;

public class Main {
    public static void main(String[] args) throws Exception {
        // Set switches to be used
        CommandLineParser parser = new DefaultParser();
        Options options = new Options();

        options.addOption("d", "descriptors", false, "Calculate descriptors");
        options.addOption("3D", false, "Compute 3D descriptors");
        options.addOption("f", "fingerprint", true, "Calculate a specific fingerprint. " +
                "Must be one of {FP, ExtFP, EStateFP, GraphFP, MACCSFP, PubchemFP, SubFP, " +
                "KRFP, AP2DFP, SubFPC, KRFPC, AP2DFPC}");
        options.addOption("nBits", true, "Number of bits of FP and GraphFP fingerprints (default: 1024)");
        options.addOption("searchDepth", true, "Search depth of FP and GraphFP fingerprints (default: 7)");
        options.addOption("n", "names", false, "Obtain only names of descriptors/fingerprint bits");
        options.addOption("i", "input", true, "Input v2000 SD file (ignored if --names)");
        //options.addOption("o", "output", false, "Output tab-separated file (ignored if --names)");
        options.addOption("h", "help", false, "Shows this Help");

        try {
            // Parse switches
            CommandLine commandLine = parser.parse(options, args);
            // Define descriptors
            List<CDK_Descriptor> descriptors = new ArrayList<>(Arrays.asList(new CDK_AcidicGroupCountDescriptor(),
                    new CDK_ALOGPDescriptor(), new CDK_APolDescriptor(), new CDK_AromaticAtomsCountDescriptor(),
                    new CDK_AromaticBondsCountDescriptor(), new CDK_AtomCountDescriptor(new String[] { "*" }),
                    new CDK_HeavyAtomCountDescriptor(), new CDK_AtomCountDescriptor(new String[] { "H" }),
                    new CDK_AtomCountDescriptor(new String[] { "B" }), new CDK_AtomCountDescriptor(new String[] { "C" }),
                    new CDK_AtomCountDescriptor(new String[] { "N" }), new CDK_AtomCountDescriptor(new String[] { "O" }),
                    new CDK_AtomCountDescriptor(new String[] { "S" }), new CDK_AtomCountDescriptor(new String[] { "P" }),
                    new CDK_HalogenCountDescriptor(), new CDK_AutocorrelationDescriptor(), new CDK_BaryszMatrixDescriptor(),
                    new CDK_BasicGroupCountDescriptor(), new CDK_BCUTDescriptor(), new CDK_BondCountDescriptor(),
                    new CDK_BPolDescriptor(), new CDK_BurdenModifiedEigenvaluesDescriptor(),
                    new CDK_CarbonTypesDescriptor(), new CDK_ChiChainDescriptor(), new CDK_ChiClusterDescriptor(),
                    new CDK_ChiPathClusterDescriptor(), new CDK_ChiPathDescriptor(), new CDK_ConstitutionalDescriptor(),
                    new CDK_CrippenDescriptor(), new CDK_DetourMatrixDescriptor(),
                    new CDK_EccentricConnectivityIndexDescriptor(), new CDK_EStateAtomTypeDescriptor(),
                    new CDK_ExtendedTopochemicalAtomDescriptor(), new CDK_FMFDescriptor(),
                    new CDK_FragmentComplexityDescriptor(), new CDK_HBondAcceptorCountDescriptor(),
                    new CDK_HBondDonorCountDescriptor(), new CDK_HybridizationRatioDescriptor(),
                    new CDK_InformationContentDescriptor(), new CDK_KappaShapeIndicesDescriptor(),
                    new CDK_LargestChainDescriptor(), new CDK_LargestPiSystemDescriptor(),
                    new CDK_LongestAliphaticChainDescriptor(), new CDK_MannholdLogPDescriptor(),
                    new CDK_McGowanVolumeDescriptor(), new CDK_MDEDescriptor(), new CDK_MLFERDescriptor(),
                    new CDK_PathCountDescriptor(), new CDK_PetitjeanNumberDescriptor(), new CDK_RingCountDescriptor(),
                    new CDK_RotatableBondsCountDescriptor(), new CDK_RuleOfFiveDescriptor(),new CDK_TopologicalDescriptor(),
                    new CDK_TopologicalChargeDescriptor(), new CDK_TopologicalDistanceMatrixDescriptor(),
                    new CDK_TPSADescriptor(), new CDK_VABCDescriptor(), new CDK_VAdjMaDescriptor(),
                    new CDK_WalkCountDescriptor(), new CDK_WeightDescriptor(), new CDK_WeightedPathDescriptor(),
                    new CDK_WienerNumbersDescriptor(), new CDK_XLogPDescriptor(), new CDK_ZagrebIndexDescriptor()));
            // Add 3D descriptors
            if (commandLine.hasOption("3D")) {
                descriptors.addAll(Arrays.asList(new CDK_Autocorrelation3DDescriptor(),
                        new CDK_CPSADescriptor(), new CDK_GravitationalIndexDescriptor(),
                        new CDK_LengthOverBreadthDescriptor(), new CDK_MomentOfInertiaDescriptor(),
                        new CDK_PetitjeanShapeIndexDescriptor(), new CDK_RDFDescriptor(), new CDK_WHIMDescriptor()));
            }
            // Fingerprint names, to be chosen from
            // Cannot instantiate as size and search depth may vary
            List<String> fingerprints = new ArrayList<>(Arrays.asList("FP", "ExtFP", "EStateFP", "GraphFP", "MACCSFP",
                    "PubchemFP", "SubFP", "KRFP", "AP2DFP", "SubFPC", "KRFPC", "AP2DFPC"));
            // Display help
            if (commandLine.hasOption("help")) {
                new HelpFormatter().printHelp("java -jar ePaDEL.jar", options);
            } else if (commandLine.hasOption("descriptors")) {
                if (commandLine.hasOption("names")) {
                    // Output names of descriptors
                    List<String> names = new ArrayList<String>();
                    for (CDK_Descriptor desc : descriptors) {
                        names.addAll(Arrays.asList(desc.getDescriptorNames()));
                    }
                    String out = String.join(" ", names);
                    System.out.print(out);
                } else {
                    // Calculate descriptor values
                    if (!commandLine.hasOption("input")){
                        // No input given
                        throw new Exception("Input V3000 SD file must be provided.");
                    }
                    // Read content of file
                    try (FileInputStream fis = new FileInputStream(commandLine.getOptionValue("input"))) {
                        IteratingMDLReader supplier = new IteratingMDLReader(fis, DefaultChemObjectBuilder.getInstance());
                        // Iterate over molecules
                        while (supplier.hasNext()) {
                            IAtomContainer molecule = supplier.next();
                            List<String> desc_values = new ArrayList<>();
                            // Iterate over descriptors
                            for (CDK_Descriptor desc : descriptors) {
                                // Set molecule in descriptor calculator and run
                                desc.setMolecule(molecule);
                                try {
                                    desc.run();
                                    desc_values.addAll(Arrays.asList(desc.getDescriptorValues()));
                                } catch (ArrayIndexOutOfBoundsException | NullPointerException e) {
                                    // Handle exceptions
                                    String[] empty = new String[desc.getDescriptorNames().length];
                                    Arrays.fill(empty, "NaN");
                                    desc_values.addAll(Arrays.asList(empty));
                                }
                            }
                            // Print values to sdtout
                            System.out.println(String.join(" ", desc_values));
                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            } else if (commandLine.hasOption("fingerprint")) {
                // Obtain type of fingerprint
                String fp_value = commandLine.getOptionValue("fingerprint");
                if (!fingerprints.contains(fp_value)){
                    // Not supported
                    throw new Exception("Fingerprint type " + fp_value + " is not available.");
                }
                int size;
                // Get size of fingerprint
                if (commandLine.hasOption("nBits")) {
                    size = Integer.parseInt(commandLine.getOptionValue("nBits"));
                }
                else {
                    size = 1024;
                }
                // Same with search depth
                int searchDepth;
                if (commandLine.hasOption("searchDepth")) {
                    searchDepth = Integer.parseInt(commandLine.getOptionValue("searchDepth"));
                }
                else {
                    searchDepth = 7;
                }
                // Instantiate fingerprint calculator
                eCDK_IFingerprint fp;
                if (fp_value.equals("FP")) {
                    fp = new eCDK_Fingerprinter(size, searchDepth);
                } else if (fp_value.equals("ExtFP")) {
                    // Size and searchDepth not supported with CDK < 2.5
                    fp = new eCDK_ExtendedFingerprinter(); //size, searchDepth);
                } else if (fp_value.equals("EStateFP")) {
                    fp = new eCDK_EStateFingerprinter();
                } else if (fp_value.equals("GraphFP")) {
                    fp = new eCDK_GraphOnlyFingerprinter(size, searchDepth);
                } else if (fp_value.equals("MACCSFP")) {
                    fp = new eCDK_MACCSFingerprinter();
                } else if (fp_value.equals("PubchemFP")) {
                    fp = new eCDK_PubchemFingerprinter();
                } else if (fp_value.equals("SubFP")) {
                    fp = new eCDK_SubstructureFingerprinter();
                } else if (fp_value.equals("KRFP")) {
                    fp = new eCDK_KlekotaRothFingerprinter();
                } else if (fp_value.equals("AP2DFP")) {
                    fp = new eCDK_AtomPairs2DFingerprinter();
                } else if (fp_value.equals("SubFPC")) {
                    fp = new eCDK_SubstructureFingerprintCount();
                } else if (fp_value.equals("KRFPC")) {
                    fp = new eCDK_KlekotaRothFingerprintCount();
                } else if (fp_value.equals("AP2DFPC")) {
                    fp = new eCDK_AtomPairs2DFingerprintCount();
                } else {
                    // Default to eCDKFingerprinter
                    fp = new eCDK_Fingerprinter(size, searchDepth);
                }
                // Obtain FP names
                if (commandLine.hasOption("names")) {
                    // Use ethane as default molecule
                    fp.setMolecule(new SmilesParser(DefaultChemObjectBuilder.getInstance()).parseSmiles("CC"));
                    List<String> names = new ArrayList<>(Arrays.asList(fp.getDescriptorNames()));
                    String out = String.join(" ", names);
                    System.out.print(out);
                } else {
                    // Calculate fingerprint values
                    if (!commandLine.hasOption("input")){
                        // No input given
                        throw new Exception("Input V3000 SD file must be provided.");
                    }
                    // Read content of file
                    try (FileInputStream fis = new FileInputStream(commandLine.getOptionValue("input"))) {
                        IteratingMDLReader supplier = new IteratingMDLReader(fis, DefaultChemObjectBuilder.getInstance());
                        // Iterate over molecules
                        while (supplier.hasNext()) {
                            IAtomContainer molecule = supplier.next();
                            List<String> fp_values = new ArrayList<>();
                            // Set molecule in descriptor calculator and run
                            fp.setMolecule(molecule);
                            try {
                                fp.run();
                                fp_values.addAll(Arrays.asList(fp.getDescriptorValues()));
                            } catch (ArrayIndexOutOfBoundsException e) {
                                // Handle exceptions
                                String[] empty = new String[fp.getDescriptorNames().length];
                                Arrays.fill(empty, "NaN");
                                fp_values.addAll(Arrays.asList(empty));
                            }
                            // Print values to sdtout
                            System.out.println(String.join(" ", fp_values));
                        }
                    } catch (Exception e) {}
                }
            }
        }
        catch (ParseException e) {
            e.printStackTrace();
        }
    }
}