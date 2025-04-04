/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_decomp.h
 * @ingroup CONSHDLRS-GCG
 * @brief  constraint handler for structure detection
 * @author Martin Bergner
 * @author Michael Bastubbe
 * @author Hanna Franzen
 * @author Erik Muehmer
 *
 * This constraint handler manages the structure detection process. It will run all registered structure detectors in an
 * iterative refinement scheme. Afterwards some post-processing detectors might be called.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CONS_DECOMP_H__
#define GCG_CONS_DECOMP_H__

#include "scip/scip.h"
#include "scip/type_retcode.h"

#include "gcg/gcg.h"


#ifdef __cplusplus
extern "C" {
#endif

/*
 * incomplete type used in C code
 */
typedef struct PARTIALDECOMP_C PARTIALDECOMP_C;

/** @brief Gets the character of the detector
 * @returns detector character */
GCG_EXPORT
char GCGdetectorGetChar(
   GCG_DETECTOR*         detector            /**< pointer to detector */
   );

/**
 * @brief returns the data of the provided detector
 * @returns data of the provided detector
 */
GCG_EXPORT
GCG_DETECTORDATA* GCGdetectorGetData(
   GCG_DETECTOR* detector  /**< Detector data structure */
   );

/**
 * @brief returns the name of the provided detector
 * @returns name of the given detector
 */
GCG_EXPORT
const char* GCGdetectorGetName(
   GCG_DETECTOR* detector  /**< detector data structure */
   );

/** @brief interface method to detect the structure including presolving
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGdetectStructure(
   GCG*                  gcg,               /**< GCG data structure */
   SCIP_RESULT*          result             /**< Result pointer to indicate whether some structure was found */
   );

/**
 * @brief searches for the consclassifier with the given name and returns it or NULL if classifier is not found
 * @returns consclassifier pointer or NULL if consclassifier with given name is not found
 */
GCG_EXPORT
GCG_CONSCLASSIFIER* GCGfindConsClassifier(
   GCG* gcg,            /**< GCG data structure  */
   const char* name     /**< the name of the searched consclassifier */
   );

/**
 * @brief searches for the varclassifier with the given name and returns it or NULL if classifier is not found
 * @returns varclassifier pointer or NULL if varclassifier with given name is not found
 */
GCG_EXPORT
GCG_VARCLASSIFIER* GCGfindVarClassifier(
   GCG* gcg,            /**< GCG data structure  */
   const char* name     /**< the name of the searched varclassifier */
   );


/**
 * @brief searches for the detector with the given name and returns it or NULL if detector is not found
 * @returns detector pointer or NULL if detector with given name is not found
 */
GCG_EXPORT
GCG_DETECTOR* GCGfindDetector(
   GCG* gcg,            /**< GCG data structure  */
   const char* name     /**< the name of the searched detector */
   );

/** @brief method to adapt score for orig decomps
 * @todo: change score for some parameter settings
 *
 * @returns new score */
GCG_EXPORT
SCIP_Real GCGconshdlrDecompAdaptScore(
   GCG*              gcg,     /**< GCG data structure */
   SCIP_Real         oldscore /**< current score (to be updated) */
   );

/**
 * @brief searches for the score with the given name and returns it or NULL if score is not found
 * @returns score pointer or NULL if score with given name is not found
 */
GCG_SCORE* GCGconshdlrDecompFindScore(
   GCG*                  gcg,
   const char*           name
   );

/**
 * @brief searches for the score with the given shortname and returns it or NULL if score is not found
 * @returns score pointer or NULL if score with given shortname is not found
 */
GCG_SCORE* GCGconshdlrDecompFindScoreByShortname(
   GCG*                  gcg,
   const char*           shortname
   );

/** @brief Gets the best known decomposition
 *
 * @note caller has to free returned GCG_DECOMP
 * @returns the decomposition if available and NULL otherwise */
GCG_EXPORT
GCG_DECOMP* GCGgetBestDecomp(
   GCG*                  gcg,                /**< GCG data structure */
   SCIP_Bool             printwarnings       /**< should warnings pre printed */
   );

/**
 * @brief returns the remaining time of scip that the decomposition may use
 * @returns remaining  time that the decompositon may use
 */
GCG_EXPORT
SCIP_Real GCGgetRemainingTime(
   GCG*                 gcg                 /**< GCG data structure */
   );

/**
 * @brief includes one constraint classifier
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGincludeConsClassifier(
   GCG*                  gcg,             /**< GCG data structure */
   const char*           name,            /**< name of the classifier */
   const char*           description,     /**< describing main idea of this classifier */
   int                   priority,        /**< priority of the classifier */
   SCIP_Bool             enabled,         /**< whether the classifier should be enabled by default */
   GCG_CLASSIFIERDATA*   classifierdata,  /**< classifierdata the associated classifier data (or NULL) */
   GCG_DECL_FREECONSCLASSIFIER((*freeClassifier)),  /**< destructor of classifier (or NULL) */
   GCG_DECL_CONSCLASSIFY((*classify))               /**< the method that will classify constraints or variables (must not be NULL) */
   );

/**
 * @brief includes one detector
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGincludeDetector(
   GCG*                  gcg,                     /**< GCG data structure */
   const char*           name,                    /**< name of the detector */
   const char            decchar,                 /**< char that is used in detector chain history for this detector */
   const char*           description,             /**< describing main idea of this detector */
   int                   freqCallRound,           /**< frequency the detector gets called in detection loop, i.e. it is called in round r if and only if minCallRound <= r <= maxCallRound AND (r - minCallRound) mod freqCallRound == 0 */
   int                   maxCallRound,            /**< last detection round the detector gets called */
   int                   minCallRound,            /**< first round the detector gets called (offset in detection loop) */
   int                   freqCallRoundOriginal,   /**< frequency the detector gets called in detection loop while detecting of the original problem */
   int                   maxCallRoundOriginal,    /**< last round the detector gets called while detecting of the original problem */
   int                   minCallRoundOriginal,    /**< first round the detector gets called (offset in detection loop) while detecting of the original problem */
   int                   priority,                /**< priority of the detector */
   SCIP_Bool             enabled,                 /**< whether the detector should be enabled by default */
   SCIP_Bool             enabledFinishing,        /**< whether the finishing should be enabled */
   SCIP_Bool             enabledPostprocessing,   /**< whether the postprocessing should be enabled */
   SCIP_Bool             skip,                    /**< whether the detector should be skipped if others found structure */
   SCIP_Bool             usefulRecall,            /**< is it useful to call this detector on a descendant of the propagated partialdec */
   GCG_DETECTORDATA      *detectordata,           /**< the associated detector data (or NULL) */
   GCG_DECL_FREEDETECTOR((*freeDetector)),        /**< destructor of detector (or NULL) */
   GCG_DECL_INITDETECTOR((*initDetector)),        /**< initialization method of detector (or NULL) */
   GCG_DECL_EXITDETECTOR((*exitDetector)),        /**< deinitialization method of detector (or NULL) */
   GCG_DECL_PROPAGATEPARTIALDEC((*propagatePartialdecDetector)),      /**< method to refine a partial decomposition inside detection loop (or NULL) */
   GCG_DECL_FINISHPARTIALDEC((*finishPartialdecDetector)),            /**< method to complete a partial decomposition when called in detection loop (or NULL) */
   GCG_DECL_POSTPROCESSPARTIALDEC((*postprocessPartialdecDetector)),  /**< method to postprocess a complete decomposition, called after detection loop (or NULL) */
   GCG_DECL_SETPARAMAGGRESSIVE((*setParamAggressiveDetector)),        /**< method that is called if the detection emphasis setting aggressive is chosen */
   GCG_DECL_SETPARAMDEFAULT((*setParamDefaultDetector)),              /**< method that is called if the detection emphasis setting default is chosen */
   GCG_DECL_SETPARAMFAST((*setParamFastDetector))                     /**< method that is called if the detection emphasis setting fast is chosen */
   );

/**
 * @brief includes one variable classifier
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGincludeVarClassifier(
   GCG*                  gcg,           /**< GCG data structure */
   const char*           name,          /**< name of the classifier */
   const char*           description,   /**< description of the classifier */
   int                   priority,      /**< priority how early classifier is invoked */
   SCIP_Bool             enabled,       /**< whether the classifier should be enabled by default */
   GCG_CLASSIFIERDATA*   classifierdata,/**< classifierdata the associated classifier data (or NULL) */
   GCG_DECL_FREEVARCLASSIFIER((*freeClassifier)),   /**< destructor of classifier (or NULL) */
   GCG_DECL_VARCLASSIFY((*classify))                /**< method that will classify variables (must not be NULL) */
   );

/**
 * @brief includes one score
 * @returns scip return code
 */
SCIP_RETCODE GCGconshdlrDecompIncludeScore(
   GCG*                  gcg,                /**< GCG data structure */
   const char*           name,               /**< name of the score */
   const char*           shortname,          /**< shortname of the score */
   const char*           description,        /**< description of the score */
   GCG_SCOREDATA*        scoredata,          /**< scoredata the associated score data (or NULL) */
   GCG_DECL_SCOREFREE    ((*scorefree)),     /**< destructor of score (or NULL) */
   GCG_DECL_SCORECALC    ((*scorecalc))      /**< method that will calculate the scorevalue (must not be NULL) */
   );

/** @brief writes out a list of all detectors */
GCG_EXPORT
void GCGprintListOfDetectors(
   GCG*                  gcg                 /**< GCG data structure */
   );

/** @brief write out all known decompositions
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGwriteAllDecomps(
   GCG*                  gcg,                /**< GCG data structure */
   char*                 directory,          /**< directory for decompositions */
   char*                 extension,          /**< the file extension for the export */
   SCIP_Bool             original,           /**< should decomps for original problem be written */
   SCIP_Bool             presolved           /**< should decomps for preoslved problem be written */
   );

/** @brief writes all selected decompositions
 * @returns scip return code
*/
GCG_EXPORT
SCIP_RETCODE GCGwriteSelectedDecomps(
   GCG*                  gcg,                /**< GCG data structure */
   char*                 directory,          /**< directory for decompositions */
   char*                 extension           /**< extension for decompositions */
   );

/** @brief adds a candidate for block number and counts how often a candidate is added */
GCG_EXPORT
void GCGconshdlrDecompAddCandidatesNBlocks(
   GCG* gcg,                    /**< GCG data structure */
   SCIP_Bool origprob,           /**< which DETPROBDATA that should be modified */
   int candidate                 /**< proposed amount of blocks */
);

/**
 * @brief creates and adds a basic partialdecomp (all cons/vars are assigned to master)
 *
 * @returns SCIP_RETCODE
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompAddBasicPartialdec(
   GCG* gcg,                        /**< GCG data structure */
   SCIP_Bool presolved,             /**< create basic partialdecomp for presolved if true, otherwise for original */
   PARTIALDECOMP_C** partialdecomp  /**< pointer to the C pointer of the created partialdecomp */
   );

/**
 * @brief creates a pure matrix partialdecomp (i.e. all cons/vars to one single block)
 *
 * matrix is added to list of all partialdecs
 * @returns id of matrix partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompAddMatrixPartialdec(
   GCG* gcg,            /**< GCG data structure */
   SCIP_Bool presolved  /**< create matrix for presolved if true, otherwise for original */
   );

/**
 * @brief adds a decomp that exists before the detection is called
 * @note this method should only be called if there is no partialdec for this decomposition
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompAddPreexistingDecomp(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** @brief adds a candidate for block size given by the user */
GCG_EXPORT
void GCGconshdlrDecompAddUserCandidatesNBlocks(
   GCG* gcg,                    /**< GCG data structure */
   int candidate                 /**< candidate for block size */
   );

/**
 * @brief checks if two pricing problems are identical based on information from detection
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompArePricingprobsIdenticalForPartialdec(
   GCG*                    gcg,                 /**< GCG data structure */
   PARTIALDECOMP_C*        partialdecomp,       /**< partial decompostion for which the pricing problems are checked for identity */
   int                     probnr1,             /**< index of first block to check */
   int                     probnr2,             /**< index of second block to check */
   SCIP_Bool*              identical            /**< bool pointer to score the result of the check*/
   );

/**
 * @brief calculates and adds block size candidates using constraint classifications and variable classifications
 */
GCG_EXPORT
void GCGconshdlrDecompCalcCandidatesNBlocks(
   GCG* gcg,               /**< GCG data structure */
   SCIP_Bool transformed   /**< whether to find the candidates for the transformed problem, otherwise the original */
);

/**
 * @brief check whether partialdecs are consistent
 *
 * Checks whether
 *  1) the predecessors of all finished partialdecs in both detprobdatas can be found
 *  2) selected list is synchron with selected information in partialdecs
 *  3) selected exists is synchronized with selected list
 *
 *  @returns true if partialdec information is consistent */
 GCG_EXPORT
SCIP_Bool GCGconshdlrDecompCheckConsistency(
   GCG* gcg   /**< GCG data structure **/
   );

/**
 * @brief run classification of vars and cons
 *
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompClassify(
   GCG*                 gcg,           /**< GCG data structure */
   SCIP_Bool            transformed    /**< whether to classify the transformed problem, otherwise the original */
);

/**
 * @brief for two identical pricing problems a corresponding varmap is created
 * @param scip scip data structure
 * @param hashorig2pricingvar mapping from orig to pricingvar
 * @param partialdecomp partial decompostion for which the pricing problems are checked for identity
 * @param probnr1 index of first block
 * @param probnr2 index of second block
 * @param scip1 subscip of first block
 * @param scip2 subscip of second block
 * @param varmap mapping from orig to pricingvar
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompCreateVarmapForPartialdec(
   GCG*                    gcg,
   SCIP_HASHMAP**          hashorig2pricingvar,
   PARTIALDECOMP_C*        partialdecomp,
   int                     probnr1,
   int                     probnr2,
   SCIP*                   scip1,
   SCIP*                   scip2,
   SCIP_HASHMAP*           varmap
   );

/**
 * @brief decreases the counter for created decompositions and returns it
 * @returns number of created decompositions that was recently decreased
 */
GCG_EXPORT
int GCGconshdlrDecompDecreaseNCallsCreateDecomp(
   GCG*                 gcg                 /**< GCG data structure **/
   );

/** @brief deregisters partialdecs in the conshdlr
 *
 * Use this function for deletion of ALL the partialdecs.
 */
GCG_EXPORT
void GCGconshdlrDecompDeregisterPartialdecs(
   GCG* gcg,   /**< GCG data structure */
   SCIP_Bool original  /**< iff TRUE the status with respect to the original problem is returned */
   );

/** @brief Frees Detprobdata of the original and transformed/presolved problem.
 *
 * @note Does not free Detprobdata of the original problem if GCGconshdlrDecompFreeOrigOnExit is set to false.
 */
void GCGconshdlrDecompFreeDetprobdata(
   GCG* gcg   /**< GCG data structure */
   );

/**
 * @brief sets freeing of detection data of original problem during exit to true
 *
 * used before calling SCIPfreeTransform(),
 * set to true to revoke presolving
 * (e.g. if unpresolved decomposition is used, and transformation is not successful)
 */
GCG_EXPORT
void GCGconshdlrDecompFreeOrigOnExit(
   GCG* gcg,     /**< GCG data structure */
   SCIP_Bool free /**< whether to free orig data */
   );

/**
 * @brief returns block number user candidate with given index
 * @param scip SCIP data structure
 * @param index index of block number user candidate that should be returned
 * @returns block number user candidate with given index
 */
GCG_EXPORT
 int GCGconshdlrDecompGetBlockNumberCandidate(
    GCG*                  gcg,
    int                   index
     );

/**
 * @brief returns the total detection time
 * @param scip SCIP data structure
 * @returns total detection time
 */
GCG_EXPORT
SCIP_Real GCGconshdlrDecompGetCompleteDetectionTime(
    GCG*                  gcg
    );

/** @brief returns an array containing all decompositions
 *
 *  Updates the decomp decomposition structure by converting all finished partialdecs into decompositions and replacing the
 *  old list in the conshdlr.
 *
 *  @returns decomposition array
 *   */
GCG_EXPORT
GCG_DECOMP** GCGconshdlrDecompGetDecomps(
   GCG* gcg   /**< GCG data structure */
   );

/** @brief Gets an array of all detectors
 *
 * @returns array of detectors */
GCG_EXPORT
GCG_DETECTOR** GCGconshdlrDecompGetDetectors(
   GCG* gcg   /**< GCG data structure */
   );

/** @brief Gets an array of all scores
 *
 * @returns array of scores */
GCG_SCORE** GCGconshdlrDecompGetScores(
   GCG* gcg
   );

/**
 * @brief Gets the shortname of the currently enabled score
 * @returns the shortname of the currently enabled score
 */
GCG_EXPORT
char* GCGgetCurrentScoreShortname(
   GCG*                 gcg
   );

/**
 * @brief Gets the currently enabled score
 * @returns the currently enabled score
 */
GCG_EXPORT
GCG_SCORE* GCGgetCurrentScore(
   GCG*                 gcg
   );

/** @brief gets an array of all constraint classifier
 *
 * @returns array of constraint classifier */
GCG_EXPORT
GCG_CONSCLASSIFIER** GCGconshdlrDecompGetConsClassifiers(
   GCG* gcg   /**< GCG data structure */
   );

/** @brief gets an array of all variable classifier
 *
 * @returns array of variable classifier */
GCG_EXPORT
GCG_VARCLASSIFIER** GCGconshdlrDecompGetVarClassifiers(
   GCG* gcg
   );

/** @brief Gets a list of ids of the current partialdecs that are finished
 *
 *  @note recommendation: when in doubt plan for as many ids as partialdecs
 *  @see GCGconshdlrDecompGetNPartialdecs
 *  @returns scip return code */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompGetFinishedPartialdecsList(
   GCG*           gcg,        /**< GCG data structure */
   int**          idlist,     /**< id list to output to */
   int*           listlength  /**< length of output list */
   );

/** @brief Gets a list of ids of the current partialdecs
 *
 *  @note recommendation: when in doubt plan for as many ids as partialdecs
 *  @see GCGconshdlrDecompGetNPartialdecs
 *  @returns scip return code */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompGetPartialdecsList(
   GCG*           gcg,        /**< GCG data structure */
   int**          idlist,     /**< id list to output to */
   int*           listlength  /**< length of output list */
);

/**
 * @brief returns the number of block candidates given by the user
 * @returns number of block candidates given by the user
 */
GCG_EXPORT
 int GCGconshdlrDecompGetNBlockNumberCandidates(
   GCG*                 gcg                 /**< GCG data structure */
    );

/** @brief gets block number of partialdec
 * @returns block number of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompPartialdecGetNBlocks(
   PARTIALDECOMP_C*        partialdecomp       /**< partialdec */
   );

/** @brief gets the number of decompositions (= amount of finished partialdecs)
 *
 * @returns number of decompositions */
GCG_EXPORT
int GCGconshdlrDecompGetNDecomps(
   GCG* gcg   /**< GCG data structure */
   );

/** @brief Gets the number of all detectors
 * @returns number of detectors */
GCG_EXPORT
int GCGconshdlrDecompGetNDetectors(
   GCG* gcg   /**< GCG data structure */
   );

/** @brief Gets the number of all scores
 * @returns number of scores */
int GCGconshdlrDecompGetNScores(
   GCG* gcg   /**< GCG data structure */
   );

/** @brief Gets the number of all constraint classifiers
 * @returns number of constraint classifiers */
GCG_EXPORT
int GCGconshdlrDecompGetNConsClassifiers(
   GCG* gcg   /**< GCG data structure */
   );

/** @brief Gets the number of all variable classifiers
 * @returns number of variable classifiers */
GCG_EXPORT
int GCGconshdlrDecompGetNVarClassifiers(
   GCG* gcg
   );

/** @brief Gets the next partialdec id managed by cons_decomp
 * @returns the next partialdec id managed by cons_decomp */
GCG_EXPORT
int GCGconshdlrDecompGetNextPartialdecID(
   GCG*   gcg    /**< GCG data structure **/
   );

/** @brief gets number of active constraints during the detection of the decomp with given id
 *
 * Gets the number of constraints that were active while detecting the decomposition originating from the partialdec with the
 * given id, this method is used to decide if the problem has changed since detection, if so the aggregation information
 * needs to be recalculated
 *
 * @note if the partialdec is not complete the function returns -1
 *
 * @returns number of constraints that were active while detecting the decomposition
 */
GCG_EXPORT
int GCGconshdlrDecompGetNFormerDetectionConssForID(
   GCG* gcg,     /**< GCG data structure */
   int id         /**< id of the partialdec the information is asked for */
   );

/** @brief returns whether aggregation information has been calculated for a partialdec
  * @returns bool that indicates whether aggregation information has been calculated
 */
GCG_EXPORT
SCIP_Bool GCGconshdlrDecompPartialdecAggregationInformationCalculated(
   PARTIALDECOMP_C*        partialdecomp        /**< partialdec */
   );

/** @brief calculates aggregation information for a partialdec
 */
GCG_EXPORT
void GCGconshdlrDecompPartialdecCalcAggregationInformation(
   PARTIALDECOMP_C*        partialdecomp,           /**< partialdec */
   SCIP_Bool               ignoreDetectionLimits    /**< Set to true if computation should ignore detection limits. This parameter is ignored if the patched bliss version is not present. */
   );

/** @brief gets number of equivalence classes of partialdec
 * @returns number of equivalence classes of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompPartialdecGetNEquivalenceClasses(
   PARTIALDECOMP_C*        partialdecomp        /**< partialdec */
   );

/** @brief gets the representative block of a equivalence class of partialdec
 * @returns representative block
 */
GCG_EXPORT
int GCGconshdlrDecompPartialdecGetReprBlockForEqClass(
   PARTIALDECOMP_C*        partialdecomp,       /**< partialdec */
   int                     eqclass              /**< id of the equivalence class */
   );

/** @brief gets the blocks of a equivalence class of partialdec
 * @returns representative block
 */
GCG_EXPORT
const int* GCGconshdlrDecompPartialdecGetBlocksForEqClass(
   PARTIALDECOMP_C*        partialdecomp,       /**< partialdec */
   int                     eqclass              /**< id of the equivalence class */
   );

/** @brief gets the number of blocks of a equivalence class of partialdec
 * @returns representative block
 */
GCG_EXPORT
int GCGconshdlrDecompPartialdecGetNBlocksForEqClass(
   PARTIALDECOMP_C*        partialdecomp,       /**< partialdec */
   int                     eqclass              /**< id of the equivalence class */
   );

/** @brief Returns a vector that maps probvar indices of a block contained in an equivalence class to the probvar indices of the representative block of the class
 * @returns representative block
 */
GCG_EXPORT
const int* GCGconshdlrDecompPartialdecGetRepVarMap(
   PARTIALDECOMP_C*        partialdecomp,       /**< partialdec */
   int                     eqclass,             /**< id of the equivalence class */
   int                     eqclassblock         /** index of block (with respect to eqclass) */
   );

/** @brief gets number of linking variables of partialdec
 * @returns number of linking variables of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompPartialdecGetNLinkingVars(
   PARTIALDECOMP_C*        partialdecomp        /**< partialdec */
   );

/** @brief gets number of master constraints of partialdec
 * @returns number of master constraints of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompPartialdecGetNMasterConss(
   PARTIALDECOMP_C*        partialdecomp        /**< partialdec */
   );

/** @brief gets number of master variables of partialdec
 * @returns number of master variables of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompPartialdecGetNMasterVars(
   PARTIALDECOMP_C*        partialdecomp        /**< partialdec */
   );

/** @brief gets number of open constraints of partialdec
 * @returns total number of open constraints of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompPartialdecGetNOpenConss(
   PARTIALDECOMP_C*        partialdecomp        /**< partialdec */
   );

/** @brief gets number of open variables of partialdec
 * @returns total number of open variables of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompPartialdecGetNOpenVars(
   PARTIALDECOMP_C*        partialdecomp        /**< partialdec */
   );

/** @brief returns the original variable for a block and index
 * @returns corresponding original variable
 */
GCG_EXPORT
SCIP_VAR* GCGconshdlrDecompPartialdecGetOrigVarForBlock(
   PARTIALDECOMP_C*        partialdecomp,       /**< partialdec */
   int                     block,               /**< id of the block */
   int                     blockvarindex        /**< index of var in block */
   );

/** @brief gets the number of variables of a block
 * @returns representative block
 */
GCG_EXPORT
int GCGconshdlrDecompPartialdecGetNVarsForBlock(
   PARTIALDECOMP_C*        partialdecomp,       /**< partialdec */
   int                     block                /**< id of the block */
);

/** @brief Gets the number of finished partialdecs available for the original problem
 * @returns number of partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNFinishedPartialdecsOrig(
   GCG*       gcg   /**< GCG data structure */
   );

/** @brief Gets the number of finished partialdecs available for the transformed problem
 * @returns number of partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNFinishedPartialdecsTransformed(
   GCG*       gcg   /**< GCG data structure */
   );

/** @brief Gets the number of open partialdecs available for the original problem
 * @returns number of partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNOpenPartialdecsOrig(
   GCG*       gcg   /**< GCG data structure */
);

/** @brief Gets the number of open partialdecs available for the transformed problem
 * @returns number of partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNOpenPartialdecsTransformed(
   GCG*       gcg   /**< GCG data structure */
);

/** @brief Gets the number of all partialdecs
 * @returns number of Partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNPartialdecs(
   GCG*       gcg   /**< GCG data structure */
   );

/** @brief Gets the number of partialdecs available for the original problem
 * @returns number of partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNPartialdecsOrig(
   GCG*       gcg   /**< GCG data structure */
   );

/** @brief Gets the number of partialdecs available for the transformed problem
 * @returns number of partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNPartialdecsTransformed(
   GCG*       gcg   /**< GCG data structure */
   );

/** @brief gets number of stairlinking variables of partialdec
 * @returns total number of stairlinking variables of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompPartialdecGetNStairlinkingVars(
   PARTIALDECOMP_C*        partialdecomp        /**< partialdec */
   );

/** @brief Gets wrapped PARTIALDECOMP with given id
 *
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompGetPartialdecFromID(
   GCG*                gcg,              /**< GCG data structure */
   int                  partialdecid,     /**< id of PARTIALDECOMP */
   PARTIALDECOMP_C**    partialdecomp     /**< pointer to C wrapper pointer for PARTIALDECOMP */
   );

/** @brief gets score of partialdec
 * @returns score in respect to current score type
 */
GCG_EXPORT
SCIP_Real GCGconshdlrDecompGetPartialdecScore(
   PARTIALDECOMP_C*        partialdecomp        /**< partialdec */
   );

/** @brief gets the clock tracking the score computation time
 * @returns pointer to a clock object
 */
GCG_EXPORT
SCIP_CLOCK* GCGconshdlrDecompGetScoreClock(
   GCG* gcg      /**< GCG data structure */
   );

/** @brief gets total score computation time
 * @returns total score computation time
 */
GCG_EXPORT
SCIP_Real GCGconshdlrDecompGetScoreTotalTime(
   GCG* gcg      /**< GCG data structure */
   );

/** @brief Gets a list of ids of all currently selected partialdecs
 *  @returns list of partialdecs */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompGetSelectedPartialdecs(
   GCG*          gcg,        /**< GCG data structure */
   int**          idlist,     /**< id list to output to */
   int*           listlength  /**< length of output list */
   );

/**
 * @brief counts up the counter for created decompositions and returns it
 * @returns number of created decompositions that was recently increased
 */
GCG_EXPORT
int GCGconshdlrDecompIncreaseNCallsCreateDecomp(
  GCG*                 gcg                 /**< GCG data structure **/
   );

/** @brief gets whether partialdec is presolved
 * @returns true iff partialdec is presolved
 */
GCG_EXPORT
SCIP_Bool GCGconshdlrDecompPartialdecIsPresolved(
   PARTIALDECOMP_C*        partialdecomp        /**< partialdec */
   );

/** @brief gets whether partialdec is selected
 * @returns true iff partialdec is selected
 */
GCG_EXPORT
SCIP_Bool GCGconshdlrDecompPartialdecIsSelected(
   PARTIALDECOMP_C*        partialdececomp      /**< partialdec */
   );

/**
 * @brief returns whether or not a detprobdata structure for the original problem exists
 * @returns true iff an original detprobdata exists
 */
GCG_EXPORT
SCIP_Bool GCGconshdlrDecompOrigDetprobdataExists(
   GCG*                 gcg                 /**< GCG data structure */
   );

/**
 * @brief returns whether or not an original decompositions exists in the data structures
 * @returns true iff an origial decomposition exist
 */
GCG_EXPORT
SCIP_Bool GCGconshdlrDecompOrigPartialdecExists(
   GCG*                 gcg                 /**< GCG data structure */
   );

/**
 * @brief returns whether or not a detprobdata structure for the presolved problem exists
 * @returns true iff a presolved detprobdata exists
 */
GCG_EXPORT
SCIP_Bool GCGconshdlrDecompPresolvedDetprobdataExists(
   GCG*                 gcg                 /**< GCG data structure */
   );

/** @brief display statistics about detectors
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompPrintDetectorStatistics(
   GCG*                 gcg,                /**< GCG data structure */
   FILE*                 file                /**< output file or NULL for standard output */
   );

/** @brief display statistics about scores
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompPrintScoreStatistics(
   GCG*                 gcg,                /**< GCG data structure */
   FILE*                 file                /**< output file or NULL for standard output */
   );

/**
 * @brief selects/unselects a partialdecomp
 *
 * @returns SCIP return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompSelectPartialdec(
   PARTIALDECOMP_C*        partialdecomp,       /**< partialdecomp */
   SCIP_Bool               select               /**< select/unselect */
   );

/** @brief sets detector parameters values
 *
 * sets detector parameters values to
 *
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all detector parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for detection is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the detectors produce more decompositions
 *  - SCIP_PARAMSETTING_OFF which turns off all detection
 *
 * @returns SCIP return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompSetDetection(
   GCG*                 gcg,                /**< GCG data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/**
 * @brief translates n best unpresolved partialdec to a complete presolved one
 * @param scip SCIP data structure
 * @param n number of partialdecs that should be translated
 * @param completeGreedily whether or not to complete the decomposition greedily
 * @param translateSymmetry whether or not to translate symmetry information
 * @returns SCIP return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompTranslateNBestOrigPartialdecs(
   GCG*                  gcg,
   int                   n,
   SCIP_Bool             completeGreedily,
   SCIP_Bool             translateSymmetry
);

/**
 * @brief translates unpresolved partialdec to a complete presolved one
 * @param gcg GCG data structure
 * @returns SCIP return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompTranslateOrigPartialdecs(
   GCG*                  gcg
   );

/** Gets whether the detection already took place
 * @returns true if detection took place, false otherwise */
GCG_EXPORT
SCIP_Bool GCGdetectionTookPlace(
   GCG*  gcg, /**< GCG data structure */
   SCIP_Bool original /**< iff TRUE the status with respect to the original problem is returned */
   );

/**
 * method to eliminate duplicate constraint names and name unnamed constraints
 * @return SCIP return code
 */
GCG_EXPORT
SCIP_RETCODE SCIPconshdlrDecompRepairConsNames(
   GCG*                gcg   /**< GCG data structure */
   );

/**
 * @brief creates the constraint handler for decomp and includes it in SCIP
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGincludeConshdlrDecomp(
   GCG*  gcg   /**< GCG data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
