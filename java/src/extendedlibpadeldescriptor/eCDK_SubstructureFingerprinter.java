
package extendedlibpadeldescriptor;

import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.ExtendedFingerprinter;
import libpadeldescriptor.CDK_SubstructureFingerprinter;

public class eCDK_SubstructureFingerprinter extends CDK_SubstructureFingerprinter
{
    public eCDK_SubstructureFingerprinter() {}
    
    public String[] getDescriptorValues() {
        return this.descriptorValues_.clone();
    }
}
