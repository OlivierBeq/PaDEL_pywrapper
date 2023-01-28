
This folder contains additional files required to get access to fingerprints' values.

## Compile

The `javac` and `jar` executables are assumed to be in a folder declared in your path environment variable.

Compile the additional files like so (tested with Java Development Kit 1.6.0): 

```cmd
mkdir build

javac -d build -classpath "../src/PaDEL_pywrapper/PaDEL-Descriptor/lib/libPaDEL-Descriptor.jar" src/extendedlibpadeldescriptor/*.java

jar cvf elibPaDEL-Descriptor.jar -C build .
```

Then make the compiled library available to the Python wrapper:

```cmd
cp elibPaDEL-Descriptor.jar ../src/PaDEL_pywrapper/PaDEL-Descriptor/lib/
```