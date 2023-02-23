
package extendedlibpadeldescriptor;

import org.openscience.cdk.interfaces.IAtomContainer;

public interface eCDK_IFingerprint
{
    public String[] getDescriptorValues();

    public String[] getDescriptorNames();

    public void setMolecule(IAtomContainer molecule);

    public void run();
}