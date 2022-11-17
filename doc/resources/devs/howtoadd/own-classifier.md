# How to add classifiers {#own-classifier}

[TOC]

Classifiers are called to **sort variables and constraints into groups** to use this knowledge
inside a [detector](#detectors) in a later stage of the [detection](#detection-process).
\n
A complete list of all classifiers contained in this release can be found [here](#classifiers).

# Adding your own Classifier

With the following steps, we explain how you can **add your own constraint/variable classifier plugin**:
1. **Preparations**
  1. Choose a name `myclassifier` for your classifier.
  2. Copy the template files `src/clscons_xyz.cpp`/`src/clsvar_xyz.cpp` and `src/clscons_xyz.h`/`src/clsvar_xyz.h`
   while renaming `xyz` to `myclassifier`.
  3. Open the new files with a text editor and replace all occurrences of `Xyz` by `Myclassifier` and `xyz` by `myclassifier`.
2. **Creating your Classifier**
  1. Adjust the properties of the classifier (see @ref CLS_PROPERTIES).
  2. [optional] Define the classifier data (see @ref CLS_DATA).
  3. Implement the interface methods (see @ref CLS_INTERFACE).
  4. Implement the fundamental callback methods (see @ref CLS_FUNDAMENTALCALLBACKS).
  5. [optional] Implement the additional callback methods (see @ref CLS_ADDITIONALCALLBACKS).
3. **Make GCG use it**
  1. Add it to gcgplugins.c by adding
    1. the line <tt>\#include clscons_myclassifer.h</tt> / <tt>\#include clsvar_myclassifer.h</tt> in the `/* classifiers */` section.
    2. the line `SCIP_CALL( SCIPincludeConsClassifierMyclassifier(scip) );` in  the `/* Classifiers */` section.
  2. Add it to your build system:
    1. _Using Makefile:_ Add your classifier `.o` (`clscons_myclassifier.o`/`clsvar_myclassifier.o`) to the list below `LIBOBJ =` in the file `Makefile` in the root folder.
    2. _Using CMake:_ In `src/CMakeLists.txt`, add your `clscons_myclassifier.cpp`/`clsvar_myclassifier.cpp` below `set(gcgsources` and your
   `clscons_myclassifier.h`/`clsvar_myclassifier.h` below `set(gcgheaders`.


## Properties of a Classifier {#CLS_PROPERTIES}

At the top of the new file `clscons_myclassifier.cpp`/`clsvar_myclassifier.cpp`, you can find the classifier properties.
These are given as compiler defines.
The properties you have to set have the following meaning:

@par DEC_CLASSIFIERNAME: the name of classifier
This name is used in the interactive shell to address the classifier. Names have to be unique: no two classifiers may have the same name.

@par DEC_DESC: short description of classification
This string is printed as description of the classifier in the interactive shell.

@par DEC_PRIORITY: priority of classifier
At the start of the detection process, the classifiers are called in a predefined order, which is given by the priorities of those. The classifiers are called in the order of decreasing priority.

\par DEC_ENABLED: Flag to indicate whether the classifier should be enabled by default.
Disabled classifiers are not started.

@par DEC_ENABLEDORIG: classify on original problem
Set this flag to true if the classifier should classify on the original (non-presolved) problem.

@par DEC_ENABLEDPRESOLVED: classify on presolved problem
Set this flag to true if the classifier should classify on the presolved problem.

## Classifier Data {#CLS_DATA}
Defining classifier data is optional.

## Interface Methods {#CLS_INTERFACE}
At the bottom of `clscons_myclassifier.cpp`/`clsvar_myclassifier.cpp`, you can find the interface method `SCIPincludeConsClassifierXyz()`/`SCIPincludeVarClassifierXyz()`,
which also appears in `clscons_myclassifier.h`/`clsvar_myclassifier.h`.
\n
This method has to be adjusted only slightly.
It is responsible for notifying GCG (and especially cons_decomp.c) of the presence of the classifier by calling the method
`GCGincludeConsClassifier()`/`GCGincludeVarClassifier()`.
`SCIPincludeConsClassifierXyz()`/`SCIPincludeVarClassifierXyz()` is called by the user to include the classifier,
i.e., to use the classifier in the application (see 3.1.1. at the top of the page).

If you are using classifier data, you have to allocate the memory for the data at this point.
You can do this by calling
```C
SCIP_CALL( SCIPallocMemory(scip, &classifierdata) );
```
For freeing the classifier data, see @ref CLS_FREE.

You may also add user parameters for your classifier, see the parameters documentation of SCIP for how to add user parameters.


## Fundamental Callback Methods of a Classifier {#CLS_FUNDAMENTALCALLBACKS}
The fundamental callback methods of the plug-ins are the ones that have to be implemented in order to obtain
an operational algorithm. Classifier plug-ins have one main function:
 * @ref CLS_CONSCLASSIFY "classify constraints" according to some property  
 * @ref CLS_VARCLASSIFY "classify variables" according to some property

Exactly one of following methods has to be implemented for every classifier.

Additional documentation for the callback methods of classifiers can be found in the
files type_consclassifier.h and type_varclassifier.h.

### CONSCLASSIFY {#CLS_CONSCLASSIFY}
The `DEC_DECL_CONSCLASSIFY(classifierClassify)` callback assigns constraints to classes using the `assignConsToClass()` method of the `gcg::ConsClassifier`.

### VARCLASSIFY {#CLS_VARCLASSIFY}
The `DEC_DECL_VARCLASSIFY(classifierClassify)` callback  assigns variables to classes using the `assignVarToClass()` method of the `gcg::VarClassifier`.

## Additional Callback Methods of a Classifier {#CLS_ADDITIONALCALLBACKS}
### FREECLASSIFIER {#CLS_FREE}
The `DEC_DECL_FREECLASSIFIER(classifierFreeXyz)` callback is called upon exiting GCG to free user data.
