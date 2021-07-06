# How to add decomposition scores {#own-score}

GCG uses scores to **rate decompositions** that are found by the detection. To continuously prefer different
decompositions than GCG does by default with its set partitioning foreseeing max white score, you can program your
own score.\n
A complete list of all scores contained in this release can be found @ref detection-scores "here".

# Adding your own Decomposition Score
With the following steps, we explain how you can **add your own score to rate decompositions**:
1. **Preparations**\n
  1. Choose a name `myscore` for your decomposition score.
  2. In the file class_partialdecomp.h, initialize your score:
    - In the class `PARTIALDECOMP::PARTIALDECOMP(SCIP* _scip, ...)`, initialize it below the line `bool originalProblem` by adding `myscore( -1. ),`.
    - In the class `PARTIALDECOMP::PARTIALDECOMP(const PARTIALDECOMP *partialdectocopy)`, add it to the partial decomposition to copy by adding `myscore = partialdectocopy->myscore;`
  3. [optional] In the function void PARTIALDECOMP::displayInfo(int detailLevel), add an output for the statistics for your score.
  4. Add a real variable `myscore` for your score to the `class PARTIALDECOMP` (below the section ` /* score values */`).

2. **Creating your Score**\n
In the file cons_decomp.cpp:
  1. Write your score calculation function with the following signature (as an example, you can have a look at the `GCGconshdlrDecompCalcMaxWhiteScore()` function):
  ```
  SCIP_RETCODE GCGconshdlrDecompCalcMyScore(
   SCIP* scip,
   int partialdecid,
   SCIP_Real* score
   )
  ``` 
  2. In the case switch in the function `PARTIALDECOMP::getScore(type)` in the file class_partialdecomp.h, add your calculation function.
  3. Add a getter and a setter function `getMyScore()` and `setMyScore()` (as an example, take the function `getMaxWhiteScore()`) in the file class_partialdecomp.h.

3. **Make GCG use it**\n
Add your score to 
  1. the arrays `scoretype_shortnames` and `scoretype_descriptions` in the file scoretype.c.
  2. the enum `scoretype` in the file type_scoretype.h.
  3. the `struct Dec_Scores` in the file pub_decomp.h.