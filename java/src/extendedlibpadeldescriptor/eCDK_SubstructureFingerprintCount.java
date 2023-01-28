
package extendedlibpadeldescriptor;

import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import libpadeldescriptor.CDK_SubstructureFingerprintCount;

public class eCDK_SubstructureFingerprintCount extends CDK_SubstructureFingerprintCount
{
    public eCDK_SubstructureFingerprintCount() {}
    
    public String[] getDescriptorValues() {
        return this.descriptorValues_.clone();
    }
}
