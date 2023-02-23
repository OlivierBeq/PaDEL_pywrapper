
package extendedlibpadeldescriptor;

import org.openscience.cdk.interfaces.IAtomContainer;
import libpadeldescriptor.CDK_KlekotaRothFingerprinter;

public class eCDK_KlekotaRothFingerprinter extends CDK_KlekotaRothFingerprinter implements eCDK_IFingerprint
{
    public eCDK_KlekotaRothFingerprinter() {}

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
