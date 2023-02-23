
package extendedlibpadeldescriptor;

import org.openscience.cdk.interfaces.IAtomContainer;
import libpadeldescriptor.CDK_AtomPairs2DFingerprintCount;

public class eCDK_AtomPairs2DFingerprintCount extends CDK_AtomPairs2DFingerprintCount implements eCDK_IFingerprint
{
    public eCDK_AtomPairs2DFingerprintCount() {}

    public String[] getDescriptorValues() {
        return this.descriptorValues_.clone();
    }

    public String[] getDescriptorNames() {
        String[] names = this.cdkDescriptor_.getDescriptorNames();
        for (int i=0; i < names.length; i++){
            names[i] = names[i].replace("APC2D", "AP2DFPC");
        }
        return names;
    }

    public void setMolecule(IAtomContainer molecule){
        super.setMolecule(molecule);
    }

    public void run(){
        super.run();
    }
}
