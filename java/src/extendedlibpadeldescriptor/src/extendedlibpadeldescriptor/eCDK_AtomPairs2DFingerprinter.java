
package extendedlibpadeldescriptor;

import org.openscience.cdk.interfaces.IAtomContainer;
import libpadeldescriptor.CDK_AtomPairs2DFingerprinter;

public class eCDK_AtomPairs2DFingerprinter extends CDK_AtomPairs2DFingerprinter implements eCDK_IFingerprint
{
    public eCDK_AtomPairs2DFingerprinter() {
        this.setPrefix("AP2DFP");
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
