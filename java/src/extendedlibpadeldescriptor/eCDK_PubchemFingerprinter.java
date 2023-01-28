
package extendedlibpadeldescriptor;

import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import libpadeldescriptor.CDK_PubchemFingerprinter;

public class eCDK_PubchemFingerprinter extends CDK_PubchemFingerprinter
{
    public eCDK_PubchemFingerprinter() {}
    
    public String[] getDescriptorValues() {
        return this.descriptorValues_.clone();
    }
}
