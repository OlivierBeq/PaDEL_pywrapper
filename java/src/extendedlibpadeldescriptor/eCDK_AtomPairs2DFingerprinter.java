
package extendedlibpadeldescriptor;

import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import libpadeldescriptor.CDK_AtomPairs2DFingerprinter;

public class eCDK_AtomPairs2DFingerprinter extends CDK_AtomPairs2DFingerprinter
{
    public eCDK_AtomPairs2DFingerprinter() {}
    
    public String[] getDescriptorValues() {
        return this.descriptorValues_.clone();
    }
}
