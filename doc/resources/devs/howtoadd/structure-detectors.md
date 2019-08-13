# Detectors (deprecated) {#structure-detectors}

Structure detectors are used to detect or enforce a structure suitable for Dantzig-Wolfe Reformulation (DWR).
\n
A complete list of all detectors/enforcers contained in this release can be found "here" (please add).

In the following, we explain how the user can add its own structure enforcement plug-in.
Take the basic detector (dec_connected.c) as an example.
As all other default plug-ins, it is written in C. There is currently no C++ wrapper available.

Additional documentation for the callback methods of structure detectors, in particular for their input parameters,
can be found in the file type_detector.h.

Here is what you have to do to implement a detector:
- Copy the template files src/dec_xyz.c and src/dec_xyz.h into files named "dec_mydetector.c"
   and "dec_mydetector.h".
   \n
   Make sure to adjust your Makefile such that these files are compiled and linked to your project.
- Open the new files with a text editor and replace all occurrences of "xyz" by "mydetector".
- Adjust the properties of the detector (see \ref DEC_PROPERTIES).
- Define the detector data (see \ref DEC_DATA). This is optional.
- Implement the interface methods (see \ref DEC_INTERFACE).
- Implement the fundamental callback methods (see ???).
- Implement the additional callback methods (see \ref DEC_ADDITIONALCALLBACKS). This is optional.


# Properties of a Detector {#DEC_PROPERTIES}

At the top of the new file "dec_mydetector.c", you can find the detector properties.
These are given as compiler defines.
The properties you have to set have the following meaning:

\par DEC_DETECTORNAME: the name of the detector.
This name is used in the interactive shell to address the detector.
Additionally, if you are searching for a detector with SCIPfindDetector(), this name is looked up.
Names have to be unique: no two detectors may have the same name.

\par DEC_DESC: the description of the detector.
This string is printed as description of the detector in the interactive shell.

\par DEC_PRIORITY: the priority of the detector.
At the begin of the solving process, the detectors are called in a predefined order, which is given by the priorities
of the detectors.
The detectors are called in the order of decreasing priority.
\n
The priority of the detector should be set according to the complexity of the detection algorithm and the quality of the decomposition:
detectors that provide fast algorithms that usually have a good decomposition (i.e., provide good dual bound) should have a high
priority. An easy way to list the priorities of all detectors is "display detectors" in the interactive shell of GCG.

\par DEC_DECCHAR: Display character of the detector.
The unique display character for the detector. It can be used to quickly distinguish similar structures by the detector which found them.

\par DEC_ENABLED: Flag to indicate whether the detector should be enabled by default.
Disabled detectors are not started.

\par DEC_SKIP: Flag to indicate whether the detector should be skipped if other detectors found decompositions
This flag is useful if the detector acts as a last resort to generate a decomposition. It will not be called if any detector of higher
priority found a decomposition.

# Detector Data {#DEC_DATA}

Below the header "Data structures" you can find the struct "struct DEC_DetectorData".
In this data structure, you can store the data of your detector. For example, you should store the adjustable parameters
of the detector in this data structure.
\n
Defining detector data is optional. You can leave this struct empty.


# Interface Methods {#DEC_INTERFACE}

At the bottom of "dec_mydetector.c", you can find the interface method SCIPincludeDetectionMyDetector(),
which also appears in "dec_mydetector.h".
\n
This method has to be adjusted only slightly.
It is responsible for notifying GCG (and especially cons_decomp.c) of the presence of the detector by calling the method
DECincludeDetector().
SCIPincludeDetectionMyDetector() is called by the user, if he wants to include the detector,
i.e., if he wants to use the detector in his application.

If you are using detector data, you have to allocate the memory for the data at this point.
You can do this by calling
```C
SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
```
You also have to initialize the fields in struct SCIP_DetectorData afterwards. For freeing the
detector data, see ???.

You may also add user parameters for your detector, see the parameters documentation of \SCIP for how to add user parameters and
the method SCIPincludeDetectionBorderheur() in dec_connected.c for an example.


# Fundamental Callback Methods of a Detector {#DEC_FUNDAMENTALCALLBACKS'}

The fundamental callback methods of the plug-ins are the ones that have to be implemented in order to obtain
an operational algorithm. Detector plug-ins have only one fundamental callback method, namely the DETECTSTRUCTURE method.
This method has to be implemented for every detector; the other callback methods are optional.

Additional documentation to the callback methods, in particular to their input parameters,
can be found in type_detector.h.

## DETECTSTRUCTURE

The DETECTSTRUCTURE callback is called during the detection loop and should perform the actual detection.
It should inspect the problem instance at hand and deduct some structure from the constraint matrix.
It needs to store the structure information in ??? and needs to allocate the array where to store the
information.

Typical methods called by a detector are, for example, SCIPgetVars(), SCIPGetConss(), DECcreateDecompFromMasterconss(), etc. .

# Additional Callback Methods of a Detector {#DEC_ADDITIONALCALLBACKS}

## DETECTORINIT

The INITDETECTOR callback is executed after the problem was transformed.
The detector may, e.g., use this call to initialize his detector data.
The difference between the original and the transformed problem is explained in
"What is this thing with the original and the transformed problem about?" on ???.

## DETECTOREXIT

If you are using detection data (see \ref DEC_DATA and \ref DEC_INTERFACE), you have to implement this method in order to free the detection data.
This can be done by the following procedure:
```C
static
DEC_DECL_EXITDETECTOR(decExitMydetector)
{
   DEC_DETECTORDATA* detectordata;

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   SCIPfreeMemory(scip, &detectordata);

   DECdetectorSetData(detector, NULL);

   return SCIP_OKAY;
}
```
If you have allocated memory for fields in your detector data, remember to free this memory
before freeing the detector data itself.
The DETECTOREXIT callback is executed before the solution process is started.
In this method, the detector should free all resources that have been allocated for the detection process in ???.
