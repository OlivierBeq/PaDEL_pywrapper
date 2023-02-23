
package extendedlibpadeldescriptor;

import org.openscience.cdk.interfaces.IAtomContainer;
import libpadeldescriptor.CDK_SubstructureFingerprintCount;

public class eCDK_SubstructureFingerprintCount extends CDK_SubstructureFingerprintCount implements eCDK_IFingerprint
{
    public eCDK_SubstructureFingerprintCount() {}

    public String[] getDescriptorValues() {
        return this.descriptorValues_.clone();
    }

    public String[] getDescriptorNames() {
        return this.cdkDescriptor_.getDescriptorNames();
    }

    public void setMolecule(IAtomContainer molecule){
        super.setMolecule(molecule);
    }

    public void run(){
        super.run();
    }
}
