
package extendedlibpadeldescriptor;

import org.openscience.cdk.interfaces.IAtomContainer;
import libpadeldescriptor.CDK_EStateFingerprinter;

public class eCDK_EStateFingerprinter extends CDK_EStateFingerprinter implements eCDK_IFingerprint
{
    public eCDK_EStateFingerprinter() {}

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
