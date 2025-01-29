/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_decomp.h
 * @ingroup CONSHDLRS-GCG
 * @brief  constraint handler for structure detection
 * @author Martin Bergner
 * @author Michael Bastubbe
 * @author Hanna Franzen
 *
 * This constraint handler manages the structure detection process. It will run all registered structure detectors in an
 * iterative refinement scheme. Afterwards some post-processing detectors might be called.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CONS_DECOMP_H__
#define GCG_CONS_DECOMP_H__

#include "scip/scip.h"
#include "def.h"
#include "type_detector.h"
#include "type_varclassifier.h"
#include "type_consclassifier.h"
#include "type_score.h"


#ifdef __cplusplus
extern "C" {
#endif


/** forward declarations */
struct Partialdecomp_Wrapper;
typedef struct Partialdecomp_Wrapper PARTIALDECOMP_WRAPPER;

struct Detprobdata_Wrapper;

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
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_RESULT*          result             /**< Result pointer to indicate whether some structure was found */
   );

/**
 * @brief searches for the consclassifier with the given name and returns it or NULL if classifier is not found
 * @returns consclassifier pointer or NULL if consclassifier with given name is not found
 */
GCG_EXPORT
GCG_CONSCLASSIFIER* GCGfindConsClassifier(
   SCIP* scip,          /**< SCIP data structure  */
   const char* name     /**< the name of the searched consclassifier */
   );

/**
 * @brief searches for the varclassifier with the given name and returns it or NULL if classifier is not found
 * @returns varclassifier pointer or NULL if varclassifier with given name is not found
 */
GCG_EXPORT
GCG_VARCLASSIFIER* GCGfindVarClassifier(
   SCIP* scip,          /**< SCIP data structure  */
   const char* name     /**< the name of the searched varclassifier */
   );


/**
 * @brief searches for the detector with the given name and returns it or NULL if detector is not found
 * @returns detector pointer or NULL if detector with given name is not found
 */
GCG_EXPORT
GCG_DETECTOR* GCGfindDetector(
   SCIP* scip,          /**< SCIP data structure  */
   const char* name     /**< the name of the searched detector */
   );

/**
 * @brief searches for the score with the given name and returns it or NULL if score is not found
 * @returns score pointer or NULL if score with given name is not found
 */
GCG_SCORE* GCGconshdlrDecompFindScore(
   SCIP*                 scip,
   const char*           name
   );

/**
 * @brief searches for the score with the given shortname and returns it or NULL if score is not found
 * @returns score pointer or NULL if score with given shortname is not found
 */
GCG_SCORE* GCGconshdlrDecompFindScoreByShortname(
   SCIP*                 scip,
   const char*           shortname
   );

/** @brief Gets the best known decomposition
 *
 * @note caller has to free returned GCG_DECOMP
 * @returns the decomposition if available and NULL otherwise */
GCG_EXPORT
GCG_DECOMP* GCGgetBestDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             printwarnings       /**< should warnings pre printed */
   );

/** @brief Gets the currently considered best partialdec
 *
 * If there is a partialdec marked to be returned (e.g. by /GCGwriteAllDecomps), it is written.
 * Else, the currently "best" decomp is returned.
 *
 * @returns partialdec to write if one can be found, or partialdecwrapper->partialdec = NULL otherwise */
GCG_EXPORT
SCIP_RETCODE GCGgetPartialdecToWrite(
   SCIP*                   scip,                /**< SCIP data structure */
   SCIP_Bool               transformed,         /**< is the problem transformed yet */
   PARTIALDECOMP_WRAPPER*  partialdecwrapper    /**< partialdec wrapper to output */
   );

/**
 * @brief returns the remaining time of scip that the decomposition may use
 * @returns remaining  time that the decompositon may use
 */
GCG_EXPORT
SCIP_Real GCGgetRemainingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**
 * @brief includes one constraint classifier
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGincludeConsClassifier(
   SCIP*                 scip,            /**< scip data structure */
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
   SCIP*                 scip,                    /**< scip data structure */
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
   SCIP*                 scip,          /**< scip data structure */
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
   SCIP*                 scip,               /**< scip data structure */
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
   SCIP*                 scip                /**< SCIP data structure */
   );

/** @brief write out all known decompositions
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGwriteAllDecomps(
   SCIP*                 scip,               /**< SCIP data structure */
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
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 directory,          /**< directory for decompositions */
   char*                 extension           /**< extension for decompositions */
   );

/** @brief adds a candidate for block number and counts how often a candidate is added */
GCG_EXPORT
void GCGconshdlrDecompAddCandidatesNBlocks(
   SCIP* scip,                   /**< SCIP data structure */
   SCIP_Bool origprob,           /**< which DETPROBDATA that should be modified */
   int candidate                 /**< proposed amount of blocks */
);

/**
 * @brief adds the given decomposition structure
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompAddDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp,             /**< GCG_DECOMP data structure */
   SCIP_Bool             select              /**< select the decomposition as candidate */
   );

/**
 * @brief creates and adds a basic partialdecomp (all cons/vars are assigned to master)
 *
 * @returns id of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompAddBasicPartialdec(
   SCIP* scip,          /**< SCIP data structure */
   SCIP_Bool presolved  /**< create basic partialdecomp for presolved if true, otherwise for original */
   );

/**
 * @brief creates a pure matrix partialdecomp (i.e. all cons/vars to one single block)
 *
 * matrix is added to list of all partialdecs
 * @returns id of matrix partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompAddMatrixPartialdec(
   SCIP* scip,          /**< SCIP data structure */
   SCIP_Bool presolved  /**< create matrix for presolved if true, otherwise for original */
   );

/**
 * @brief adds a decomp that exists before the detection is called
 * @note this method should only be called if there is no partialdec for this decomposition
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompAddPreexistingDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_DECOMP*           decomp              /**< decomposition data structure */
   );

/** @brief adds a candidate for block size given by the user */
GCG_EXPORT
void GCGconshdlrDecompAddUserCandidatesNBlocks(
   SCIP* scip,                   /**< SCIP data structure */
   int candidate                 /**< candidate for block size */
   );

/**
 * @brief checks if two pricing problems are identical based on information from detection
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompArePricingprobsIdenticalForPartialdecid(
   SCIP*                scip,             /**< scip scip data structure */
   int                  partialdecid,     /**< partialdecid id of the partial decompostion for which the pricing problems are checked for identity */
   int                  probnr1,          /**< index of first block to check */
   int                  probnr2,          /**< index of second block to check */
   SCIP_Bool*           identical         /**< bool pointer to score the result of the check*/
   );

/**
 * @brief calculates and adds block size candidates using constraint classifications and variable classifications
 */
GCG_EXPORT
void GCGconshdlrDecompCalcCandidatesNBlocks(
   SCIP* scip,             /**< SCIP data structure */
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
   SCIP* scip  /**< SCIP data structure **/
   );

/**
 * @brief run classification of vars and cons
 *
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompClassify(
   SCIP*                scip,          /**< SCIP data structure */
   SCIP_Bool            transformed    /**< whether to classify the transformed problem, otherwise the original */
);

/**
 * @brief for two identical pricing problems a corresponding varmap is created
 * @param scip scip data structure
 * @param hashorig2pricingvar mapping from orig to pricingvar
 * @param partialdecid id of the partial decompostion for which the pricing problems are checked for identity
 * @param probnr1 index of first block
 * @param probnr2 index of second block
 * @param scip1 subscip of first block
 * @param scip2 subscip of second block
 * @param varmap mapping from orig to pricingvar
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompCreateVarmapForPartialdecId(
   SCIP*                scip,
   SCIP_HASHMAP**       hashorig2pricingvar,
   int                  partialdecid,
   int                  probnr1,
   int                  probnr2,
   SCIP*                scip1,
   SCIP*                scip2,
   SCIP_HASHMAP*        varmap
   );

/**
 * @brief decreases the counter for created decompositions and returns it
 * @returns number of created decompositions that was recently decreased
 */
GCG_EXPORT
int GCGconshdlrDecompDecreaseNCallsCreateDecomp(
   SCIP*                 scip                /**< SCIP data structure **/
   );

/** @brief deregisters partialdecs in the conshdlr
 *
 * Use this function for deletion of ALL the partialdecs.
 */
GCG_EXPORT
void GCGconshdlrDecompDeregisterPartialdecs(
   SCIP* scip,  /**< SCIP data structure */
   SCIP_Bool original  /**< iff TRUE the status with respect to the original problem is returned */
   );

/** @brief Frees Detprobdata of the original and transformed/presolved problem.
 *
 * @note Does not free Detprobdata of the original problem if GCGconshdlrDecompFreeOrigOnExit is set to false.
 */
void GCGconshdlrDecompFreeDetprobdata(
   SCIP* scip  /**< SCIP data structure */
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
   SCIP* scip,    /**< SCIP data structure */
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
    SCIP*                 scip,
    int                   index
     );

/**
 * @brief returns the total detection time
 * @param scip SCIP data structure
 * @returns total detection time
 */
GCG_EXPORT
SCIP_Real GCGconshdlrDecompGetCompleteDetectionTime(
    SCIP*                 scip
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
   SCIP* scip  /**< SCIP data structure */
   );

/** @brief Gets an array of all detectors
 *
 * @returns array of detectors */
GCG_EXPORT
GCG_DETECTOR** GCGconshdlrDecompGetDetectors(
   SCIP* scip  /**< SCIP data structure */
   );

/** @brief Gets an array of all scores
 *
 * @returns array of scores */
GCG_SCORE** GCGconshdlrDecompGetScores(
   SCIP* scip
   );

/**
 * @brief Gets the shortname of the currently enabled score
 * @returns the shortname of the currently enabled score
 */
GCG_EXPORT
char* GCGgetCurrentScoreShortname(
   SCIP*                 scip
   );

/**
 * @brief Gets the currently enabled score
 * @returns the currently enabled score
 */
GCG_EXPORT
GCG_SCORE* GCGgetCurrentScore(
   SCIP*                 scip
   );

/** @brief gets an array of all constraint classifier
 *
 * @returns array of constraint classifier */
GCG_EXPORT
GCG_CONSCLASSIFIER** GCGconshdlrDecompGetConsClassifiers(
   SCIP* scip  /**< SCIP data structure */
   );

/** @brief gets an array of all variable classifier
 *
 * @returns array of variable classifier */
GCG_EXPORT
GCG_VARCLASSIFIER** GCGconshdlrDecompGetVarClassifiers(
   SCIP* scip
   );

/** @brief Gets a list of ids of the current partialdecs that are finished
 *
 *  @note recommendation: when in doubt plan for as many ids as partialdecs
 *  @see GCGconshdlrDecompGetNPartialdecs
 *  @returns scip return code */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompGetFinishedPartialdecsList(
   SCIP*          scip,       /**< SCIP data structure */
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
   SCIP*          scip,       /**< SCIP data structure */
   int**          idlist,     /**< id list to output to */
   int*           listlength  /**< length of output list */
);

/**
 * @brief returns the number of block candidates given by the user
 * @returns number of block candidates given by the user
 */
GCG_EXPORT
 int GCGconshdlrDecompGetNBlockNumberCandidates(
   SCIP*                 scip                /**< SCIP data structure */
    );

/** @brief gets block number of partialdec with given id
 * @returns block number of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompGetNBlocksByPartialdecId(
   SCIP* scip,    /**< SCIP data structure */
   int id         /**< id of partialdec */
   );

/** @brief gets the number of decompositions (= amount of finished partialdecs)
 *
 * @returns number of decompositions */
GCG_EXPORT
int GCGconshdlrDecompGetNDecomps(
   SCIP* scip  /**< SCIP data structure */
   );

/** @brief Gets the number of all detectors
 * @returns number of detectors */
GCG_EXPORT
int GCGconshdlrDecompGetNDetectors(
   SCIP* scip  /**< SCIP data structure */
   );

/** @brief Gets the number of all scores
 * @returns number of scores */
int GCGconshdlrDecompGetNScores(
   SCIP* scip  /**< SCIP data structure */
   );

/** @brief Gets the number of all constraint classifiers
 * @returns number of constraint classifiers */
GCG_EXPORT
int GCGconshdlrDecompGetNConsClassifiers(
   SCIP* scip  /**< SCIP data structure */
   );

/** @brief Gets the number of all variable classifiers
 * @returns number of variable classifiers */
GCG_EXPORT
int GCGconshdlrDecompGetNVarClassifiers(
   SCIP* scip
   );

/** @brief Gets the next partialdec id managed by cons_decomp
 * @returns the next partialdec id managed by cons_decomp */
GCG_EXPORT
int GCGconshdlrDecompGetNextPartialdecID(
   SCIP*   scip   /**< SCIP data structure **/
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
   SCIP* scip,    /**< SCIP data structure */
   int id         /**< id of the partialdec the information is asked for */
   );

/** @brief gets number of linking variables of partialdec with given id
 * @returns number of linking variables of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompGetNLinkingVarsByPartialdecId(
   SCIP* scip,    /**< SCIP data structure */
   int id         /**< id of partialdec */
   );

/** @brief gets number of master constraints of partialdec with given id
 * @returns number of master constraints of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompGetNMasterConssByPartialdecId(
   SCIP* scip,    /**< SCIP data structure */
   int id         /**< id of partialdec */
   );

/** @brief gets number of master variables of partialdec with given id
 * @returns number of master variables of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompGetNMasterVarsByPartialdecId(
   SCIP* scip,    /**< SCIP data structure */
   int id         /**< id of partialdec */
   );

/** @brief gets number of open constraints of partialdec with given id
 * @returns total number of open constraints of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompGetNOpenConssByPartialdecId(
   SCIP* scip,    /**< SCIP data structure */
   int id         /**< id of partialdec */
   );

/** @brief gets number of open variables of partialdec with given id
 * @returns total number of open variables of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompGetNOpenVarsByPartialdecId(
   SCIP* scip,    /**< SCIP data structure */
   int id         /**< id of partialdec */
   );

/** @brief Gets the number of finished partialdecs available for the original problem
 * @returns number of partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNFinishedPartialdecsOrig(
   SCIP*       scip  /**< SCIP data structure */
   );

/** @brief Gets the number of finished partialdecs available for the transformed problem
 * @returns number of partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNFinishedPartialdecsTransformed(
   SCIP*       scip  /**< SCIP data structure */
   );

/** @brief Gets the number of open partialdecs available for the original problem
 * @returns number of partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNOpenPartialdecsOrig(
   SCIP*       scip  /**< SCIP data structure */
);

/** @brief Gets the number of open partialdecs available for the transformed problem
 * @returns number of partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNOpenPartialdecsTransformed(
   SCIP*       scip  /**< SCIP data structure */
);

/** @brief Gets the number of all partialdecs
 * @returns number of Partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNPartialdecs(
   SCIP*       scip  /**< SCIP data structure */
   );

/** @brief Gets the number of partialdecs available for the original problem
 * @returns number of partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNPartialdecsOrig(
   SCIP*       scip  /**< SCIP data structure */
   );

/** @brief Gets the number of partialdecs available for the transformed problem
 * @returns number of partialdecs */
GCG_EXPORT
unsigned int GCGconshdlrDecompGetNPartialdecsTransformed(
   SCIP*       scip  /**< SCIP data structure */
   );

/** @brief gets number of stairlinking variables of partialdec with given id
 * @returns total number of stairlinking variables of partialdec
 */
GCG_EXPORT
int GCGconshdlrDecompGetNStairlinkingVarsByPartialdecId(
   SCIP* scip,    /**< SCIP data structure */
   int id         /**< id of partialdec */
   );

/** @brief Gets wrapped PARTIALDECOMP with given id
 *
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompGetPartialdecFromID(
   SCIP*          scip,             /**< SCIP data structure */
   int            partialdecid,     /**< id of PARTIALDECOMP */
   PARTIALDECOMP_WRAPPER* pwr       /**< wrapper for output PARTIALDECOMP */
   );

/** @brief gets score of partialdec with given id
 * @returns score in respect to current score type
 */
GCG_EXPORT
SCIP_Real GCGconshdlrDecompGetScoreByPartialdecId(
   SCIP* scip,    /**< SCIP data structure */
   int id         /**< id of partialdec */
   );

/** @brief gets the clock tracking the score computation time
 * @returns pointer to a clock object
 */
GCG_EXPORT
SCIP_CLOCK* GCGconshdlrDecompGetScoreClock(
   SCIP* scip     /**< SCIP data structure */
   );

/** @brief gets total score computation time
 * @returns total score computation time
 */
GCG_EXPORT
SCIP_Real GCGconshdlrDecompGetScoreTotalTime(
   SCIP* scip     /**< SCIP data structure */
   );

/** @brief Gets a list of ids of all currently selected partialdecs
 *  @returns list of partialdecs */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompGetSelectedPartialdecs(
   SCIP*          scip,       /**< SCIP data structure */
   int**          idlist,     /**< id list to output to */
   int*           listlength  /**< length of output list */
   );

/**
 * @brief counts up the counter for created decompositions and returns it
 * @returns number of created decompositions that was recently increased
 */
GCG_EXPORT
int GCGconshdlrDecompIncreaseNCallsCreateDecomp(
  SCIP*                 scip                /**< SCIP data structure **/
   );

/** @brief gets whether partialdec with given id is presolved
 * @returns true iff partialdec is presolved
 */
GCG_EXPORT
SCIP_Bool GCGconshdlrDecompIsPresolvedByPartialdecId(
   SCIP* scip,    /**< SCIP data structure */
   int id         /**< id of partialdec */
   );

/** @brief gets whether partialdec with given id is selected
 * @returns true iff partialdec is selected
 */
GCG_EXPORT
SCIP_Bool GCGconshdlrDecompIsSelectedByPartialdecId(
   SCIP* scip,    /**< SCIP data structure */
   int id         /**< id of partialdec */
   );

/**
 * @brief returns whether or not a detprobdata structure for the original problem exists
 * @returns true iff an original detprobdata exists
 */
GCG_EXPORT
SCIP_Bool GCGconshdlrDecompOrigDetprobdataExists(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**
 * @brief returns whether or not an original decompositions exists in the data structures
 * @returns true iff an origial decomposition exist
 */
GCG_EXPORT
SCIP_Bool GCGconshdlrDecompOrigPartialdecExists(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**
 * @brief returns whether or not a detprobdata structure for the presolved problem exists
 * @returns true iff a presolved detprobdata exists
 */
GCG_EXPORT
SCIP_Bool GCGconshdlrDecompPresolvedDetprobdataExists(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** @brief display statistics about detectors
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompPrintDetectorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
   );

/** @brief display statistics about scores
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompPrintScoreStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
   );

/**
 * @brief selects/unselects a partialdecomp
 *
 * @returns SCIP return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompSelectPartialdec(
   SCIP* scip,          /**< SCIP data structure */
   int partialdecid,    /**< id of partialdecomp */
   SCIP_Bool select     /**< select/unselect */
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
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/**
 * @brief translates n best unpresolved partialdec to a complete presolved one
 * @param scip SCIP data structure
 * @param n number of partialdecs that should be translated
 * @param completeGreedily whether or not to complete the decomposition greedily
 * @returns SCIP return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompTranslateNBestOrigPartialdecs(
   SCIP*                 scip,
   int                   n,
   SCIP_Bool             completeGreedily
);

/**
 * @brief translates unpresolved partialdec to a complete presolved one
 * @param scip SCIP data structure
 * @returns SCIP return code
 */
GCG_EXPORT
SCIP_RETCODE GCGconshdlrDecompTranslateOrigPartialdecs(
   SCIP*                 scip
   );

/** Gets whether the detection already took place
 * @returns true if detection took place, false otherwise */
GCG_EXPORT
SCIP_Bool GCGdetectionTookPlace(
   SCIP*  scip, /**< SCIP data structure */
   SCIP_Bool original /**< iff TRUE the status with respect to the original problem is returned */
   );

/**
 * method to eliminate duplicate constraint names and name unnamed constraints
 * @return SCIP return code
 */
GCG_EXPORT
SCIP_RETCODE SCIPconshdlrDecompRepairConsNames(
   SCIP*                scip  /**< SCIP data structure */
   );

/**
 * @brief creates the constraint handler for decomp and includes it in SCIP
 * @returns scip return code
 */
GCG_EXPORT
SCIP_RETCODE SCIPincludeConshdlrDecomp(
   SCIP* scip  /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
