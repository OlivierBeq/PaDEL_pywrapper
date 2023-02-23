
package extendedlibpadeldescriptor;

import org.openscience.cdk.interfaces.IAtomContainer;
import libpadeldescriptor.CDK_PubchemFingerprinter;

public class eCDK_PubchemFingerprinter extends CDK_PubchemFingerprinter implements eCDK_IFingerprint
{
    public eCDK_PubchemFingerprinter() {}

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
