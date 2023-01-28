
package extendedlibpadeldescriptor;

import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import libpadeldescriptor.CDK_AtomPairs2DFingerprintCount;

public class eCDK_AtomPairs2DFingerprintCount extends CDK_AtomPairs2DFingerprintCount
{
    public eCDK_AtomPairs2DFingerprintCount() {}
    
    public String[] getDescriptorValues() {
        return this.descriptorValues_.clone();
    }
}
