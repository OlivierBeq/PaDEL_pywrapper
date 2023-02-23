
package extendedlibpadeldescriptor;

import org.openscience.cdk.interfaces.IAtomContainer;
import libpadeldescriptor.CDK_SubstructureFingerprinter;

public class eCDK_SubstructureFingerprinter extends CDK_SubstructureFingerprinter implements eCDK_IFingerprint
{
    public eCDK_SubstructureFingerprinter() {}

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
