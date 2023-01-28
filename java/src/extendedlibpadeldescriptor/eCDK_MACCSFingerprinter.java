
package extendedlibpadeldescriptor;

import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import libpadeldescriptor.CDK_MACCSFingerprinter;

public class eCDK_MACCSFingerprinter extends CDK_MACCSFingerprinter
{
    public eCDK_MACCSFingerprinter() {}
    
    public String[] getDescriptorValues() {
        return this.descriptorValues_.clone();
    }
}
