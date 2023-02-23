
package extendedlibpadeldescriptor;

import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.fingerprint.GraphOnlyFingerprinter;
import libpadeldescriptor.CDK_Fingerprint;

public class eCDK_GraphOnlyFingerprinter extends CDK_Fingerprint implements eCDK_IFingerprint
{
    public eCDK_GraphOnlyFingerprinter(int size, int searchDepth) {
        this.cdkFingerprinter_ = (IFingerprinter)new GraphOnlyFingerprinter(size, searchDepth);
        this.initDescriptorsValues();
        this.setPrefix("GraphFP");
    }

    public eCDK_GraphOnlyFingerprinter(int size) {
        this.cdkFingerprinter_ = (IFingerprinter)new GraphOnlyFingerprinter(size);
        this.initDescriptorsValues();
        this.setPrefix("GraphFP");
    }

    public eCDK_GraphOnlyFingerprinter() {
        this.cdkFingerprinter_ = (IFingerprinter)new GraphOnlyFingerprinter();
        this.initDescriptorsValues();
        this.setPrefix("GraphFP");
    }

    public String[] getDescriptorValues() {
        return this.descriptorValues_.clone();
    }

    public String[] getDescriptorNames() {
        return super.getDescriptorNames();
    }

    public void setMolecule(IAtomContainer molecule){
        super.setMolecule(molecule);
    }

    public void run(){
        super.run();
    }
}
