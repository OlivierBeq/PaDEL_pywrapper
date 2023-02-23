
package extendedlibpadeldescriptor;

import org.openscience.cdk.interfaces.IAtomContainer;
import libpadeldescriptor.CDK_KlekotaRothFingerprintCount;

public class eCDK_KlekotaRothFingerprintCount extends CDK_KlekotaRothFingerprintCount implements eCDK_IFingerprint
{
    public eCDK_KlekotaRothFingerprintCount() {}

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
