# Your own Classifier {#own-classifier}
> **This page is currently being refactored. Some things might still be outdated.**

[TOC]

Classifiers are used to **sort variables and constraints into groups** to use this knowledge
inside a [detector](#detectors) in a later stage of the [detection](#detection-process).
\n
A complete list of all classifiers contained in this release can be found [here](#classifiers).

With the following steps, we explain how you can **add your own constraint/variable classifier plugin**:
1. Choose a name `myclassifier` for your classifier.
2. Copy the template files `src/clsons_XYZ.cpp`/`src/clsvar_XYZ.cpp` and `src/clsons_XYZ.h`/`src/clsvar_XYZ.h`
   while renaming `XYZ` to `myclassifier`.
3. Adjust your Makefile such that these files are compiled and linked to your project by adding the respective line
   to the list under `LIBOBJ =` in the file `Makefile` in the root folder.
4. Open the new files with a text editor and replace all occurrences of `XYZ` by `myclassifier`.
5. Adjust the properties of the detector (see @ref CLS_PROPERTIES).
6. [optional] Define the detector data (see @ref CLS_DATA).
7. Implement the interface methods (see @ref CLS_INTERFACE).
8. Implement the fundamental callback methods (see @ref CLS_FUNDAMENTALCALLBACKS).
9. [optional] Implement the additional callback methods (see @ref CLS_ADDITIONALCALLBACKS).


# Properties of a Classifier {#CLS_PROPERTIES}

At the top of the new file clscons_XYZ.cpp/clsvar_XYZ.cpp, you can find the detector properties.
These are given as compiler defines.
The properties you have to set have the following meaning:

@par DEC_CLASSIFIERNAME: the name of classifier
This name is used in the interactive shell to address the classifier. Names have to be unique: no two detectors may have the same name.

@par DEC_DESC: short description of classification
This string is printed as description of the classifier in the interactive shell.

@par DEC_PRIORITY: priority of classifier

@par DEC_ENABLEDORIG: classify on original problem

@par DEC_ENABLEDPRESOLVED: classify on presolved problem

# Classifier Data {#CLS_DATA}

Defining classifier data is optional.

# Interface Methods {#CLS_INTERFACE}
At the bottom of clscons_XYZ.cpp/clsvar_XYZ.cpp, you can find the interface method SCIPincludeConsClassifierXYZ()/SCIPincludeVarClassifierXYZ(),
which also appears in clscons_XYZ.h/clsvar_XYZ.h.
\n
This method has to be adjusted only slightly.
It is responsible for notifying GCG (and especially cons_decomp.c) of the presence of the detector by calling the method
DECincludeConsClassifier()/DECincludeVarClassifier().
SCIPincludeConsClassifierXYZ()/SCIPincludeVarClassifierXYZ() is called by the user to include the detector,
i.e., to use the detector in the application.

If you are using detector data, you have to allocate the memory for the data at this point.
You can do this by calling
```C
SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
```
You also have to initialize the fields in struct SCIP_DetectorData afterwards. For freeing the
detector data, see .

You may also add user parameters for your detector, see the parameters documentation of \SCIP for how to add user parameters and
the method SCIPincludeDetectionBorderheur() in dec_connected.c for an example.


# Fundamental Callback Methods of a Classifier {#CLS_FUNDAMENTALCALLBACKS}
The fundamental callback methods of the plug-ins are the ones that have to be implemented in order to obtain
an operational algorithm. Classifier plug-ins have one main function:
 * classify constraints/variables according to some property
The following method has to be implemented for every classifier.

Additional documentation for the callback methods of classifiers can be found in the
files type_consclassifier.h and type_varclassifier.h.
## DEC_DECL_CONSCLASSIFY


# Additional Callback Methods of a Classifier {#CLS_ADDITIONALCALLBACKS}

## classifierFree
## classifierInit
