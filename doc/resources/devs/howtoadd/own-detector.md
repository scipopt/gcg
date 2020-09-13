# How to add detectors {#own-detector}

[TOC]

Structure detectors are used to **detect or enforce a structure suitable for Dantzig-Wolfe Reformulation (DWR)**.
\n
A complete list of all detectors contained in this release can be found [here](#detectors).

# Adding your own Detector

With the following steps, we explain how you can **add your own structure detection plug-in**:
1. **Preparations**
  1. Choose a name `mydetector` for your detector.
  2. Copy the template files `src/dec_xyz.cpp` and `src/dec_xyz.h`
     while renaming `xyz` to `mydetector`.
  3. Open the new files with a text editor and replace all occurrences of `Xyz` by `Mydetector` and `xyz` by `mydetector`.
2. **Creating your Detector**
  1. Adjust the properties of the detector (see @ref DEC_PROPERTIES).
  2. [optional] Define the detector data (see @ref DEC_DATA).
  3. Implement the interface methods (see @ref DEC_INTERFACE).
  4. Implement the fundamental callback methods (see @ref DEC_FUNDAMENTALCALLBACKS).
  5. [optional] Implement the additional callback methods (see @ref DEC_ADDITIONALCALLBACKS).
3. **Make GCG use it**
  1. Add it to gcgplugins.c by adding
    1. the line <tt>\#include dec_mydetector.h</tt> in the `/* detection */` section.
    2. the line `SCIP_CALL( SCIPincludeDetectorMydetector(scip) );` in  the `/* Detectors and decompositions */` section.
  2. Add it to your build system:
    1. _Using Makefile:_ Adjust your Makefile such that these files are compiled and linked to your project by adding your detector with ending `.o` (`dec_mydetector.o`) to the list under `LIBOBJ =` in the file `Makefile` in the root folder.
    2. _Using CMake:_ In `src/CMakeLists.txt`, add your `dec_mydetector.cpp` below `set(gcgsources` and your `dec_mydetector.h` below the line `set(gcgheaders`.



## Properties of a Detector {#DEC_PROPERTIES}
At the top of the new file `dec_mydetector.cpp`, you can find the detector properties.
These are given as compiler defines.
The properties you have to set have the following meaning:

\par DEC_DETECTORNAME: the name of the detector.
This name is used in the interactive shell to address the detector.
Additionally, if you are searching for a detector with `SCIPfindDetector()`, this name is looked up.
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

## Detector Data {#DEC_DATA}
Below the header "Data structures" you can find the struct "struct DEC_DetectorData".
In this data structure, you can store the data of your detector. For example, you should store the adjustable parameters
of the detector in this data structure.
\n
Defining detector data is optional. You can leave this struct empty.


## Interface Methods {#DEC_INTERFACE}
At the bottom of `dec_mydetector.cpp`, you can find the interface method `SCIPincludeDetectorMydetector()`,
which also appears in `dec_mydetector.h`.
\n
This method has to be adjusted only slightly.
It is responsible for notifying GCG (and especially cons_decomp.cpp) of the presence of the detector by calling the method
`DECincludeDetector()`.
`SCIPincludeDetectorMydetector()` is called by the user to include the detector,
i.e., to use the detector in the application.

If you are using detector data, you have to allocate the memory for the data at this point.
You can do this by calling
```C
SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
```
You also have to initialize the fields in struct SCIP_DetectorData afterwards. For freeing the
detector data, see @ref DEC_ADDITIONALCALLBACKS.

You may also add user parameters for your detector, see the parameters documentation of SCIP for how to add user parameters.


## Fundamental Callback Methods of a Detector {#DEC_FUNDAMENTALCALLBACKS}
The fundamental callback methods of the plug-ins are the ones that have to be implemented in order to obtain
an operational algorithm. Detector plug-ins have three main functions:
 * Propagating (assigning variables/constraints to block or master),
 * Finishing (given a partialdec (incomplete decomposition), find finished partialdecs) and
 * Postprocessing (postprocess a given finished partialdec to find a different yet promising one).
At least one of the following methods has to be implemented for every detector; the other callback methods are optional.

Additional documentation to the callback methods, in particular to their input parameters,
can be found in type_detector.h.

### PROPAGATEPARTIALDEC
The `DEC_DECL_PROPAGATEPARTIALDEC(propagatePartialdecMydetector)` callback should assign variables/constraints to block or master. You can either create the decomposition by calling `DECcreateDecompFromMasterconss()` or `DECfilloutDecompFromConstoblock()` or use the getter and setter functions in pub_decomp.h to fill the decomposition structure.

### FINISHPARTIALDEC
The `DEC_DECL_FINISHPARTIALDEC(finishPartialdecMydetector)` callback should, given a partial decomposition, finish it.

### POSTPROCESSPARTIALDEC
The `DEC_DECL_POSTPROCESSPARTIALDEC(postprocessPartialdecMydetector)` callback should postprocess a given finished partial decomposition to find a different yet promising one.

## Additional Callback Methods of a Detector {#DEC_ADDITIONALCALLBACKS}
### INITDETECTOR {#DEC_INIT}
The `DEC_DECL_INITDETECTOR(detectorInitMydetector)` callback is executed after the problem was transformed.
The detector may, e.g., use this call to initialize his detector data.
The difference between the original and the transformed problem is explained in
@ref original-vs-transformed.

### EXITDETECTOR {#DEC_EXIT}
The `DEC_DECL_EXITDETECTOR(detectorExitMydetector)` callback has to be implemented if you are using detection data (see @ref DEC_DATA and @ref DEC_INTERFACE) in order to free the detection data.
This can be done by the following procedure:
```C
static
DEC_DECL_EXITDETECTOR(detectorExitMydetector)
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
The detectorExit callback is executed before the solution process is started.
In this method, the detector should free all resources that have been allocated for the detection process in @ref DEC_INIT.

### FREEDETECTOR {#DEC_FREE}
The destructor of the detector to free user data (called when GCG is exiting) has to be defined in `DEC_DECL_FREEDETECTOR(detectorFreeMydetector)`.

### SETPARAMAGGRESSIVE {#DEC_PARAM_AGG}
The parameters for the setting "aggressive" can be modified using the method `DEC_DECL_SETPARAMAGGRESSIVE(setParamAggressiveMydetector)`.

### SETPARAMDEFAULT {#DEC_PARAM_DEF}
The parameters for the setting "default" can be modified using the method `DEC_DECL_SETPARAMDEFAULT(setParamDefaultMydetector)`.

### SETPARAMFAST {#DEC_PARAM_FAST}
The parameters for the setting "fast" can be modified using the method `DEC_DECL_SETPARAMFAST(setParamFastMydetector)`.
