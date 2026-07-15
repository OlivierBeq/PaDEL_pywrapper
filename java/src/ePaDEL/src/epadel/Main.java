
package epadel;

import java.io.FileInputStream;
import java.io.IOException;
import java.lang.String;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.util.concurrent.TimeoutException;
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

    /**
     * Builds a fresh list of 1D/2D (and, if requested, 3D) descriptor calculators.
     * Called once per molecule when computing values, so that a molecule whose
     * calculation is abandoned after timing out (see runWithTimeout) never leaves a
     * descriptor instance shared with -- and racing against -- a later molecule.
     */
    private static List<CDK_Descriptor> createDescriptors(boolean include3D) {
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
        if (include3D) {
            descriptors.addAll(Arrays.asList(new CDK_Autocorrelation3DDescriptor(),
                    new CDK_CPSADescriptor(), new CDK_GravitationalIndexDescriptor(),
                    new CDK_LengthOverBreadthDescriptor(), new CDK_MomentOfInertiaDescriptor(),
                    new CDK_PetitjeanShapeIndexDescriptor(), new CDK_RDFDescriptor(), new CDK_WHIMDescriptor()));
        }
        return descriptors;
    }

    /**
     * Instantiates the fingerprint calculator matching fp_value. Called once per
     * molecule when computing values, for the same reason as createDescriptors above.
     */
    private static eCDK_IFingerprint createFingerprint(String fp_value, int size, int searchDepth) {
        if (fp_value.equals("FP")) {
            return new eCDK_Fingerprinter(size, searchDepth);
        } else if (fp_value.equals("ExtFP")) {
            // Size and searchDepth not supported with CDK < 2.5
            return new eCDK_ExtendedFingerprinter(); //size, searchDepth);
        } else if (fp_value.equals("EStateFP")) {
            return new eCDK_EStateFingerprinter();
        } else if (fp_value.equals("GraphFP")) {
            return new eCDK_GraphOnlyFingerprinter(size, searchDepth);
        } else if (fp_value.equals("MACCSFP")) {
            return new eCDK_MACCSFingerprinter();
        } else if (fp_value.equals("PubchemFP")) {
            return new eCDK_PubchemFingerprinter();
        } else if (fp_value.equals("SubFP")) {
            return new eCDK_SubstructureFingerprinter();
        } else if (fp_value.equals("KRFP")) {
            return new eCDK_KlekotaRothFingerprinter();
        } else if (fp_value.equals("AP2DFP")) {
            return new eCDK_AtomPairs2DFingerprinter();
        } else if (fp_value.equals("SubFPC")) {
            return new eCDK_SubstructureFingerprintCount();
        } else if (fp_value.equals("KRFPC")) {
            return new eCDK_KlekotaRothFingerprintCount();
        } else if (fp_value.equals("AP2DFPC")) {
            return new eCDK_AtomPairs2DFingerprintCount();
        } else {
            // Default to eCDKFingerprinter
            return new eCDK_Fingerprinter(size, searchDepth);
        }
    }

    /**
     * Runs a single descriptor/fingerprint calculation with a hard wall-clock timeout.
     *
     * Some PaDEL/CDK 3D descriptors (e.g. WHIM's Jama-based eigenvalue decomposition)
     * can fail to converge and loop forever on degenerate 3D geometry (very few heavy
     * atoms, near-linear/near-planar arrangements) instead of raising -- and since that
     * native numerical code never checks for interruption, it cannot be cooperatively
     * cancelled. calculation is therefore run on a throwaway daemon thread: if it has
     * not finished within timeoutMillis, a TimeoutException is thrown and the thread is
     * abandoned (it cannot be recovered, but being a daemon thread it does not stop the
     * JVM from exiting once the rest of the batch is done). Any exception the
     * calculation itself throws is re-thrown here so callers can handle it exactly as
     * if it had been run in place.
     *
     * If timeoutMillis is null, runs synchronously in the calling thread with no
     * timeout at all, matching the historical (pre-timeout) behavior.
     */
    private static void runWithTimeout(final Runnable calculation, final Long timeoutMillis) throws Exception {
        if (timeoutMillis == null) {
            calculation.run();
            return;
        }
        final Throwable[] error = new Throwable[1];
        Thread worker = new Thread(() -> {
            try {
                calculation.run();
            } catch (Throwable t) {
                error[0] = t;
            }
        });
        worker.setDaemon(true);
        worker.start();
        worker.join(timeoutMillis);
        if (worker.isAlive()) {
            throw new TimeoutException(
                    "Calculation did not complete within " + (timeoutMillis / 1000.0) + "s");
        }
        if (error[0] != null) {
            if (error[0] instanceof Exception) {
                throw (Exception) error[0];
            }
            throw new RuntimeException(error[0]);
        }
    }

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
        options.addOption("timeout", true, "Maximum number of seconds allowed for each descriptor's or " +
                "fingerprint's calculation for a single molecule; if exceeded, NaN is returned for it instead " +
                "of aborting the whole batch (default: no timeout)");
        options.addOption("h", "help", false, "Shows this Help");

        try {
            // Parse switches
            CommandLine commandLine = parser.parse(options, args);
            // Per-descriptor/fingerprint timeout, shared by both the descriptor and
            // fingerprint calculation loops below.
            Long timeoutMillis = null;
            if (commandLine.hasOption("timeout")) {
                double timeoutSeconds = Double.parseDouble(commandLine.getOptionValue("timeout"));
                timeoutMillis = Math.round(timeoutSeconds * 1000.0);
            }
            // Define descriptors
            List<CDK_Descriptor> descriptors = createDescriptors(commandLine.hasOption("3D"));
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
                        boolean include3D = commandLine.hasOption("3D");
                        // Iterate over molecules
                        while (supplier.hasNext()) {
                            IAtomContainer molecule = supplier.next();
                            List<String> desc_values = new ArrayList<>();
                            // Fresh descriptor instances per molecule: see createDescriptors.
                            for (CDK_Descriptor desc : createDescriptors(include3D)) {
                                // Set molecule in descriptor calculator and run
                                desc.setMolecule(molecule);
                                try {
                                    runWithTimeout(desc::run, timeoutMillis);
                                    desc_values.addAll(Arrays.asList(desc.getDescriptorValues()));
                                } catch (Exception e) {
                                    // Handle any per-descriptor calculation failure (e.g.
                                    // ArrayIndexOutOfBoundsException, NullPointerException,
                                    // ClassCastException from EStateAtomType on unusual atoms/elements
                                    // such as Pt in cisplatin, or a TimeoutException from
                                    // runWithTimeout) without aborting the whole batch.
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
                // Obtain FP names
                if (commandLine.hasOption("names")) {
                    // Use ethane as default molecule
                    eCDK_IFingerprint fp = createFingerprint(fp_value, size, searchDepth);
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
                            // Fresh fingerprint instance per molecule: see createFingerprint.
                            eCDK_IFingerprint fp = createFingerprint(fp_value, size, searchDepth);
                            // Set molecule in descriptor calculator and run
                            fp.setMolecule(molecule);
                            try {
                                runWithTimeout(fp::run, timeoutMillis);
                                fp_values.addAll(Arrays.asList(fp.getDescriptorValues()));
                            } catch (Exception e) {
                                // Broaden to match the descriptor loop's handling above.
                                String[] empty = new String[fp.getDescriptorNames().length];
                                Arrays.fill(empty, "NaN");
                                fp_values.addAll(Arrays.asList(empty));
                            }
                            // Print values to sdtout
                            System.out.println(String.join(" ", fp_values));
                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        }
        catch (ParseException e) {
            e.printStackTrace();
        }
    }
}
