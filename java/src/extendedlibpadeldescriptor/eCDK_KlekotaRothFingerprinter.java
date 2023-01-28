
package extendedlibpadeldescriptor;

import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import libpadeldescriptor.CDK_KlekotaRothFingerprinter;

public class eCDK_KlekotaRothFingerprinter extends CDK_KlekotaRothFingerprinter
{
    public eCDK_KlekotaRothFingerprinter() {}
    
    public String[] getDescriptorValues() {
        return this.descriptorValues_.clone();
    }
}
