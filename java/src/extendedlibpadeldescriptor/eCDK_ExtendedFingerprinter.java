
package extendedlibpadeldescriptor;

import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import libpadeldescriptor.CDK_ExtendedFingerprinter;

public class eCDK_ExtendedFingerprinter extends CDK_ExtendedFingerprinter
{
    public eCDK_ExtendedFingerprinter() {}
    
    public String[] getDescriptorValues() {
        return this.descriptorValues_.clone();
    }
}
