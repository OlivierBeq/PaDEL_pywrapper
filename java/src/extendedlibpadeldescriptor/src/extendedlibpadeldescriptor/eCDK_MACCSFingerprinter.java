
package extendedlibpadeldescriptor;

import org.openscience.cdk.interfaces.IAtomContainer;
import libpadeldescriptor.CDK_MACCSFingerprinter;

public class eCDK_MACCSFingerprinter extends CDK_MACCSFingerprinter implements eCDK_IFingerprint
{
    public eCDK_MACCSFingerprinter() {}

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
