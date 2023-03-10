
package extendedlibpadeldescriptor;

import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.fingerprint.Fingerprinter;
import libpadeldescriptor.CDK_Fingerprint;

public class eCDK_Fingerprinter extends CDK_Fingerprint implements eCDK_IFingerprint
{
    public eCDK_Fingerprinter(int size, int searchDepth) {
        this.cdkFingerprinter_ = (IFingerprinter)new Fingerprinter(size, searchDepth);
        this.initDescriptorsValues();
        this.setPrefix("FP");
    }

    public eCDK_Fingerprinter(int size) {
        this.cdkFingerprinter_ = (IFingerprinter)new Fingerprinter(size);
        this.initDescriptorsValues();
        this.setPrefix("FP");
    }

    public eCDK_Fingerprinter() {
        this.cdkFingerprinter_ = (IFingerprinter)new Fingerprinter();
        this.initDescriptorsValues();
        this.setPrefix("FP");
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
