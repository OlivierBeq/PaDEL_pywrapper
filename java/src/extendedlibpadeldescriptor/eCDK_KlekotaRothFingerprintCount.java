
package extendedlibpadeldescriptor;

import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import libpadeldescriptor.CDK_KlekotaRothFingerprintCount;

public class eCDK_KlekotaRothFingerprintCount extends CDK_KlekotaRothFingerprintCount
{
    public eCDK_KlekotaRothFingerprintCount() {}
    
    public String[] getDescriptorValues() {
        return this.descriptorValues_.clone();
    }
}
