# How to add decomposition scores {#own-score}

GCG uses scores to **rate decompositions** that are found by the detection. To continuously prefer different
decompositions than GCG does by default with its set partitioning foreseeing max white score, you can program your
own score.\n
A complete list of all scores contained in this release can be found @ref detection-scores "here".

# Adding your own Decomposition Score
With the following steps, we explain how you can **add your own score plug-in to rate decompositions**:
1. **Preparations**
  1. Choose a name `myscore` for your score.
  2. Copy the template files `src/score_xyz.cpp` and `src/score_xyz.h`
     while renaming `xyz` to `myscore`.
  3. Open the new files with a text editor and replace all occurrences of `Xyz` by `Myscore` and `xyz` by `myscore`.
2. **Creating your Score**
  1. Adjust the properties of the score (see @ref SCORE_PROPERTIES).
  2. [optional] Define the score data (see @ref SCORE_DATA).
  3. Implement the interface methods (see @ref SCORE_INTERFACE).
  4. Implement the fundamental callback methods (see @ref SCORE_FUNDAMENTALCALLBACKS).
  5. [optional] Implement the additional callback methods (see @ref SCORE_ADDITIONALCALLBACKS).
3. **Make GCG use it**
  1. Add it to gcgplugins.c by adding
    1. the line <tt>\#include score_myscore.h</tt> in the `/* scores */` section.
    2. the line `SCIP_CALL( GCGincludeScoreMyscore(scip) );` in  the `/* Scores */` section.
  2. Add it to your build system:
    1. _Using Makefile:_ Adjust your Makefile such that these files are compiled and linked to your project by adding your score with ending `.o` (`score_myscore.o`) to the list under `LIBOBJ =` in the file `Makefile` in the root folder.
    2. _Using CMake:_ In `src/CMakeLists.txt`, add your `score_myscore.cpp` below `set(gcgsources` and your `score_myscore.h` below the line `set(gcgheaders`.


## Properties of a Score {#SCORE_PROPERTIES}
At the top of the new file `score_myscore.cpp`, you can find the score properties.
These are given as compiler defines.
The properties you have to set have the following meaning:

\par SCORE_NAME: the name of the score.
This name is used in the interactive shell to address the score.
Additionally, if you are searching for a score with `GCGfindScore()`, this name is looked up.
Names have to be unique: no two scores may have the same name.

\par SCORE_SHORTNAME: the shortname of the score.
The shortname of the score is used in the explore menu to indicate the used score. 
You should use at most 6 characters for the shortname of the score.

\par SCORE_DESC: the description of the score.
This string is printed as description of the score in the interactive shell.

## Score Data {#SCORE_DATA}
Below the header "Data structures" you can find the struct "struct GCG_ScoreData".
In this data structure, you can store the data of your score. For example, you should store the adjustable parameters
of the score in this data structure.
\n
Defining score data is optional. You can leave this struct empty.

## Interface Methods {#SCORE_INTERFACE}
At the bottom of `score_myscore.cpp`, you can find the interface method `GCGincludeScoreMyscore()`,
which also appears in `score_myscore.h`.
\n
This method has to be adjusted only slightly.
It is responsible for notifying GCG (and especially cons_decomp.cpp) of the presence of the score by calling the method
`GCGincludeScore()`.
`GCGincludeScoreMyscore()` is called by the user to include the score,
i.e., to use the score in the application.

If you are using score data, you have to allocate the memory for the data at this point.
You can do this by calling
```C
SCIP_CALL( SCIPallocMemory(scip, &scoredata) );
```
You also have to initialize the fields in struct GCG_ScoreData afterwards. For freeing the
score data, see @ref SCORE_ADDITIONALCALLBACKS.

You may also add user parameters for your score, see the parameters documentation of SCIP for how to add user parameters.

## Fundamental Callback Method of a Score {#SCORE_FUNDAMENTALCALLBACKS}
The fundamental callback methods of the plug-ins are the ones that have to be implemented in order to obtain
an operational algorithm. Score plug-ins have one main function:
 * Calculating (calculating the score of a partial decomposition),
The following method has to be implemented for every score; the other callback methods are optional.

Additional documentation to the callback methods, in particular to their input parameters,
can be found in type_score.h.

### SCORECALC
The `GCG_DECL_SCORECALC(scoreCalcMyScore)` callback should calculate score for a decomposition. You can use `GCGconshdlrDecompGetPartialdecFromID` to get the a pointer of a partial decomposition.

## Additional Callback Methods of a Score {#SCORE_ADDITIONALCALLBACKS}

### SCOREFREE {#SCORE_FREE}
The destructor of the score to free user data (called when GCG is exiting) has to be defined in ` GCG_DECL_SCOREFREE(scoreFreeMyScore)`.
