
package extendedlibpadeldescriptor;

import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import libpadeldescriptor.CDK_ExtendedFingerprinter;

public class eCDK_ExtendedFingerprinter extends CDK_ExtendedFingerprinter implements eCDK_IFingerprint
{
    // Setting size and searchDepth does not work with CDK < 2.5 which is required by PaDEL
//     public eCDK_ExtendedFingerprinter(int size, int searchDepth) {
//         this.cdkFingerprinter_ = (IFingerprinter)new ExtendedFingerprinter(size, searchDepth);
//         this.initDescriptorsValues();
//         this.setPrefix("ExtFP");
//     }
//
//     public eCDK_ExtendedFingerprinter(int size) {
//         this.cdkFingerprinter_ = (IFingerprinter)new ExtendedFingerprinter(size);
//         this.initDescriptorsValues();
//         this.setPrefix("ExtFP");
//     }

    public eCDK_ExtendedFingerprinter() {
        this.cdkFingerprinter_ = (IFingerprinter)new ExtendedFingerprinter();
        this.initDescriptorsValues();
        this.setPrefix("ExtFP");
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
