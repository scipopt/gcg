/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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
 * @brief  constraint handler for structure detection
 * @author Martin Bergner
 * @author Michael Bastubbe
 *
 * This constraint handler manages the structure detection process. It will run all registered structure detectors in an
 * iterative refinement scheme. Afterwards some post-processing detectors might be called.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CONS_DECOMP_H__
#define GCG_CONS_DECOMP_H__

#include "scip/scip.h"
#include "type_detector.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief possible scores to evaluate founds decompositions
 * \sa SCIPconshdlrDecompGetScoretypeDescription for a description of this score
 * \sa SCIPconshdlrDecompGetScoretypeShortName
 * \sa class_seeed:getScore()
 */
enum scoretype {
   MAX_WHITE = 0,
   BORDER_AREA,
   CLASSIC,
   MAX_FORESSEEING_WHITE,
   SETPART_FWHITE,
   MAX_FORESEEING_AGG_WHITE,
   SETPART_AGG_FWHITE,
   BENDERS
};
typedef enum scoretype SCORETYPE;


/** forward declaration */
struct Seeed_Wrapper;
typedef struct Seeed_Wrapper SEEED_WRAPPER;


/**
 * sort the finished decompositions according to the currently chosen score in the according datastructures for the
 * presolved and original problem
 * @returns scip return code
 */
extern
SCIP_RETCODE DECconshdlrDecompSortDecompositionsByScore(
   SCIP*          scip     /**< SCIP data structure */
);


/**
 * @brief creates the constraint handler for decomp and includes it in SCIP
 * @returns scip return code
 */
extern
SCIP_RETCODE SCIPincludeConshdlrDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** \brief returns an array containing all decompositions
 *
 *  Updates the decdecomp decomposition structure by converting all finished seeeds into decompositions and replacing the
 *  old list in the conshdlr.
 *
 *  @returns decomposition array
 *   */
extern
DEC_DECOMP** SCIPconshdlrDecompGetDecdecomps(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** gets the number of decompositions (= amount of finished seeeds)
 *
 * @returns number of decompositions */
extern
int SCIPconshdlrDecompGetNDecdecomps(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**
 * Gets the number of constraints that were active while detecting the decomposition originating from the seeed with the
 * given id, this method is used to decide if the problem has changed since detection, if so the aggregation information
 * needs to be recalculated

 * @returns number of constraints that were active while detecting the decomposition
 */
extern
int SCIPconshdlrDecompGetNFormerDetectionConssForID(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   id                  /**< id of the seeed the information is asked for */
   );


/**
 * used before calling SCIPfreeTransform(),, if called to revoke presolving (e.g. if unpresolved decomposition is used, and
 * transformation is not successful), this seems mandatory to decide during consExitDecomp if the original detection
 * information should be freed
 * @returns scip return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompNotifyNonFinalFreeTransform(
   SCIP*                scip     /**< SCIP data structure */
   );


/**
 * used after calling SCIPfreeTransform() if called to revoke presolving (e.g. if unpresolved decomposition is used,
 * and transformation is not successful), this seems mandatory to decide during consExitDecomp if the original detection
 * information should be freed
 * @returns scip return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompNotifyFinishedNonFinalFreeTransform(
   SCIP*                scip     /**< SCIP data structure */
   );


/**
 * @brief returns the data of the provided detector
 * @returns data of the provided detector
 */
extern
DEC_DETECTORDATA* DECdetectorGetData(
   DEC_DETECTOR*         detector            /**< Detector data structure */
   );

/**
 * @brief returns the name of the provided detector
 * @returns name of the given detector
 */
extern
const char* DECdetectorGetName(
   DEC_DETECTOR*         detector            /**< detector data structure */
   );

/**
 * @brief searches for the detector with the given name and returns it or NULL if detector is not found
 * @returns detector pointer or NULL if detector with given name is not found
 */
extern
DEC_DETECTOR* DECfindDetector(
   SCIP*                 scip,               /**< SCIP data structure  */
   const char*           name                /**< the name of the searched detector */
   );


/**
 * @brief includes one detector
 * @param scip scip data structure
 * @param name name of the detector
 * @param decchar char that is used in detector chain history for this detector
 * @param description describing main idea of this detector
 * @param freqCallRound frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0
 * @param maxCallRound last detection round the detector gets called
 * @param minCallRound first round the detector gets called (offset in detection loop)
 * @param freqCallRoundOriginal frequency the detector gets called in detection loop while detecting of the original problem
 * @param maxCallRoundOriginal last round the detector gets called while detecting of the original problem
 * @param minCallRoundOriginal first round the detector gets called (offset in detection loop) while detecting of the original problem
 * @param priority  priority of the detector
 * @param enabled whether the detector should be enabled by default
 * @param enabledOriginal whether the detector should be enabled by default for detecting the original problem
 * @param enabledFinishing whether the finishing should be enabled
 * @param enabledPostprocessing whether the postprocessing should be enabled
 * @param skip whether the detector should be skipped if others found structure
 * @param usefulRecall is it useful to call this detector on a descendant of the propagated seeed
 * @param legacymode whether (old) DETECTSTRUCTURE method should also be used for detection
 * @param detectordata the associated detector data (or NULL)
 * @param DEC_DECL_DETECTSTRUCTURE((*detectStructure))   the method that will detect the structure (may be NULL), only used in legacy detection mode
 * @param DEC_DECL_FREEDETECTOR((*freeDetector)) destructor of detector (or NULL)
 * @param DEC_DECL_INITDETECTOR((*initDetector)) initialization method of detector (or NULL)
 * @param DEC_DECL_EXITDETECTOR((*exitDetector)) deinitialization method of detector (or NULL)
 * @param DEC_DECL_PROPAGATESEEED((*propagateSeeedDetector)) method to refine a partial decomposition inside detection loop (or NULL)
 * @param DEC_DECL_PROPAGATEFROMTOOLBOX((*propagateFromToolboxDetector)) method to refine a partial decomposition when called by user from console (or NULL)
 * @param DEC_DECL_FINISHFROMTOOLBOX((*finishFromToolboxDetector)) method to complete a partial decomposition when called by user from console (or NULL)
 * @param DEC_DECL_FINISHSEEED((*finishSeeedDetector)) method to complete a partial decomposition when called in detection loop (or NULL)
 * @param DEC_DECL_POSTPROCESSSEEED((*postprocessSeeedDetector)) method to postprocess a complete decomposition, called after detection loop (or NULL)
 * @param DEC_DECL_SETPARAMAGGRESSIVE((*setParamAggressiveDetector)) method that is called if the detection emphasis setting aggressive is chosen
 * @param DEC_DECL_SETPARAMDEFAULT((*setParamDefaultDetector))  method that is called if the detection emphasis setting default is chosen
 * @param DEC_DECL_SETPARAMFAST((*setParamFastDetector))  method that is called if the detection emphasis setting fast is chosen
 * @returns scip return code
 */
extern
SCIP_RETCODE DECincludeDetector(
   SCIP*                 scip,                             /**< SCIP data structure                                                */
   const char*           name,                             /**< name of the detector                                               */
   const char            decchar,                          /**< display character of the detector                                  */
   const char*           description,                      /**< description of the detector                                        */
   int                   freqCallRound,                    /** frequency the detector gets called in detection loop, i.e. it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
   int                   maxCallRound,                     /** last round the detector gets called                              */
   int                   minCallRound,                     /** first round the detector gets called (offset in detection loop) */
   int                   freqCallRoundOriginal,            /** frequency the detector gets called in detection loop while detecting of the original problem */
   int                   maxCallRoundOriginal,             /** last round the detector gets called while detecting of the original problem */
   int                   minCallRoundOriginal,             /** first round the detector gets called (offset in detection loop) while detecting of the original problem */
   int                   priority,                         /**< priority of the detector                                           */
   SCIP_Bool             enabled,                          /**< whether the detector should be enabled by default                  */
   SCIP_Bool             enabledOriginal,                  /**< whether the detector should be enabled by default for detecting the original problem */
   SCIP_Bool             enabledFinishing,                 /**< whether the finishing should be enabled */
   SCIP_Bool             enabledPostprocessing,            /**< whether the postprocessing should be enabled */
   SCIP_Bool             skip,                             /**< whether the detector should be skipped if others found structure   */
   SCIP_Bool             usefulRecall,                     /** is it useful to call this detector on a descendant of the propagated seeed */
   SCIP_Bool             legacymode,                       /**< whether (old) DETECTSTRUCTURE method should also be used for detection */
   DEC_DETECTORDATA      *detectordata,                    /**< the associated detector data (or NULL)                             */
   DEC_DECL_DETECTSTRUCTURE((*detectStructure)),           /**< the method that will detect the structure (must not be NULL)   */
   DEC_DECL_FREEDETECTOR((*freeDetector)),                 /**< destructor of detector (or NULL) */
   DEC_DECL_INITDETECTOR((*initDetector)),                 /**< initialization method of detector (or NULL)                        */
   DEC_DECL_EXITDETECTOR((*exitDetector)),                 /**< deinitialization method of detector (or NULL)                      */
   DEC_DECL_PROPAGATESEEED((*propagateSeeedDetector)),     /**< propagation method of detector (or NULL) */
   DEC_DECL_PROPAGATEFROMTOOLBOX((*propagateFromToolboxDetector)),   /**< propagation from toolbox method of detector (or NULL) */
   DEC_DECL_FINISHFROMTOOLBOX((*finishFromToolboxDetector)),         /**< finish from toolbox method of detector (or NULL) */
   DEC_DECL_FINISHSEEED((*finishSeeedDetector)),           /**< finish method of detector (or NULL) */
   DEC_DECL_POSTPROCESSSEEED((*postprocessSeeedDetector)), /**< postprocess method of detector (or NULL) */
   DEC_DECL_SETPARAMAGGRESSIVE((*setParamAggressiveDetector)),       /**< set method for aggressive parameters of detector (or NULL) */
   DEC_DECL_SETPARAMDEFAULT((*setParamDefaultDetector)),   /**< set method for default parameters of detector (or NULL) */
   DEC_DECL_SETPARAMFAST((*setParamFastDetector))          /**< set method for fast parameters of detector (or NULL) */
   );


/**
 * @brief returns the remaining time of scip that the decomposition may use
 * @returns remaining  time that the decompositon may use
 */
extern
SCIP_Real DECgetRemainingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );



/**
 * checks if two pricing problems are identical based on information from detection
 * @returns scip return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompArePricingprobsIdenticalForSeeedid(
   SCIP*                scip,             /**< scip scip data structure */
   int                  seeedid,          /**< seeedid id of the partial decompostion for which the pricing problems are checked for identity */
   int                  probnr1,          /**< index of first block to check */
   int                  probnr2,          /**< index of second block to check */
   SCIP_Bool*           identical         /**< bool pointer to score the result of the check*/
   );



/**
 * @brief for two identical pricing problems a corresponding varmap is created
 * @param scip scip data structure
 * @param hashorig2pricingvar  mapping from orig to pricingvar
 * @param seeedid id of the partial decompostion for which the pricing problems are checked for identity
 * @param probnr1 index of first block
 * @param probnr2 index of second block
 * @param scip1 subscip of first block
 * @param scip2 subscip of second block
 * @param varmap mapping from orig to pricingvar
 * @returns scip return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompCreateVarmapForSeeedId(
   SCIP*                scip,
   SCIP_HASHMAP**       hashorig2pricingvar,
   int                  seeedid,
   int                  probnr1,
   int                  probnr2,
   SCIP*                scip1,
   SCIP*                scip2,
   SCIP_HASHMAP*        varmap
   );


/**
 * @brief sets (and adds) the decomposition structure
 * @note this method should only be called if there is no seeed for this decomposition
 * @returns scip return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompAddDecdecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< decomposition data structure */
   );

/**
 * @brief creates the seeedpool for the presolved problem
 * @returns scip return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompCreateSeeedpool(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**
 * @brief creates the seeedpool for the unpresolved problem
 * @returns scip return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompCreateSeeedpoolUnpresolved(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**
 * @brief help method to access seeedpool for unpresolved problem
 * @returns pointer to seeed wrapper data structure
 */
extern
SEEED_WRAPPER* SCIPconshdlrDecompGetSeeedpoolUnpresolved(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**
 * @brief help method to access seeedpool for transformed problem
 * @returns pointer to seeed wrapper data structure
 */
extern
SEEED_WRAPPER* SCIPconshdlrDecompGetSeeedpool(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**
 * @brief creates a user seeed for the problem
 *  @returns SCIP return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompCreateUserSeeed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             presolved,          /**< should the user seeed be created for the presolved (transformed) problem */
   SCIP_Bool             markedincomplete    /**< boolean to notify that the created user decomposition is partial and should not be completed by assigning open constraints to the master */
   );


/**
 * @brief returns whether or not an unpresolved (untransformed) decompositions exists in the data structures
 * @returns SCIP return code
 */
extern
SCIP_Bool SCIPconshdlrDecompUnpresolvedSeeedExists(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**
 * @brief returns whether or not there exists at least one (complete or incomplete) decomposition
 * @param scip SCIP data structure
 * @returns TRUE if there exists at least one (complete or incomplete) decomposition
 */
extern
SCIP_Bool SCIPconshdlrDecompHasDecomp(
   SCIP*    scip
   );

/**
 * @brief sets the number of blocks
 *
 * set the number of blocks in the current user seeed (which is used for user input (read or modify) )
 * @returns SCIP return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetnumberOfBlocks(
   SCIP*                 scip,                /**< SCIP data structure */
   int                   nblocks              /**< number of blocks */
   );


/**
 * @brief sets a constraint by name to a block in the current user seeed
 * @param scip SCIP data structure
 * @param consname name of the constraint that should be set to a block
 * @param blockid index of the block the constraint should be assigned to
 * @returns SCIP return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsToBlock(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           consname,            /**< name of the constraint */
   int                   blockid              /**< block index ( counting from 0) */
   );

/**
 * @brief sets a constraint by name to master in the current user seeed
 * @param consname of the constraint that should be set to master
 * @returns SCIP return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsToMaster(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           consname
   );

/**
 * @brief sets a variable by name to a block in the current user seeed
 * @returns SCIP return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToBlock(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           varname,             /**< name of the variable */
   int                   blockid              /**< block index ( counting from 0) */
   );


/**
 * @brief sets a variable by name to the master in the current user seeed
 * @returns SCIP return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToMaster(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           varname              /**< name of the variable */
   );


/**
 * @brief sets a variable by name to the linking variables in the current user seeed
 * @returns SCIP return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToLinking(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           varname              /**< name of the variable */
   );

/**
 * @brief add block number user candidate (user candidates are prioritized over found ones)
 * @returns SCIP return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompAddBlockNumberCandidate(
   SCIP*                 scip,                /**< SCIP data structure */
   int                   blockNumberCandidate /**< given block number candidate */
   );


/**
 * @brief returns the number of block candidates given by the user
 * @returns number of block candidates given by the user
 */
extern
 int SCIPconshdlrDecompGetNBlockNumberCandidates(
   SCIP*                 scip                /**< SCIP data structure */
    );

 /**
  * @brief returns block number user candidate with given index
  * @param scip SCIP data structure
  * @param index index of block number user candidate that should be returned
  * @returns block number user candidate with given index
  */
extern
 int SCIPconshdlrDecompGetBlockNumberCandidate(
    SCIP*                 scip,
    int                   index
     );


 /**
  * @brief returns the total detection time
  * @param scip SCIP data structure
  * @returns total detection time
  */
 extern
 SCIP_Real SCIPconshdlrDecompGetCompleteDetectionTime(
    SCIP*                 scip
    );

 /**
  * @brief rejects and deletes the current user seeed
  * @returns SCIP return code
  */
 extern
SCIP_RETCODE SCIPconshdlrDecompUserSeeedReject(
   SCIP*                 scip                 /**< SCIP data structure */
   );


/**
 * finalizes and flushes the current user seeed, i.e. consider implicits, calc hashvalue, construct decdecomp if
 * complete etc
 * @returns SCIP return code
 */
 extern
SCIP_RETCODE SCIPconshdlrDecompUserSeeedFlush(
   SCIP*                 scip                 /**< SCIP data structure */
   );


/**
 * @brief translates unpresolved seeed to a complete presolved one
 * @param scip SCIP data structure
 * @param success  at least one unpresolved seeed could not be translated in a complete presolved one
 * @returns SCIP return code
 */
 extern
SCIP_RETCODE SCIPconshdlrDecompTranslateAndAddCompleteUnpresolvedSeeeds(
   SCIP*                 scip,
   SCIP_Bool*            success
   );


/**
 * @brief counts up the counter for created decompositions and returns it
 * @returns number of created decompositions that was recently increased
 */
 extern
int SCIPconshdlrDecompIncreaseAndGetNCallsCreateDecomp(
  SCIP*                 scip                /**< SCIP data structure **/
   );

/**
 * @brief decreases the counter for created decompositions and returns it
 * @returns number of created decompositions that was recently decreased
 */
 extern
int SCIPconshdlrDecompDecreaseAndGetNCallsCreateDecomp(
  SCIP*                 scip                /**< SCIP data structure **/
   );


/** Checks whether the currently best candidate is from the unpresolved seeedpool
 *
 * @returns true if best candidate is unpresolved, false otherwise */
 extern
SCIP_Bool SCIPconshdlrDecompIsBestCandidateUnpresolved(
   SCIP*                   scip  /**< SCIP data structure **/
   );


/** Checks whether
 *  1) the predecessors of all finished seeeds in both seeedpools can be found
 *  2) selected list is syncron with selected information in seeeds
 *  3) selected exists is synchronized with seleced list
 *
 *  @returns true if seeed information is consistent */
 extern
SCIP_Bool SCIPconshdlrDecompCheckConsistency(
   SCIP* scip  /**< SCIP data structure **/
   );

/** Gets the next seeed id managed by cons_decomp
 * @returns the next seeed id managed by cons_decomp */
 extern
int SCIPconshdlrDecompGetNextSeeedID(
   SCIP*   scip   /**< SCIP data structure **/
   );

/**
 * Gets the current scoretype
 * @returns the current scoretype */
 extern
SCORETYPE SCIPconshdlrDecompGetCurrScoretype(
   SCIP* scip  /**< SCIP data structure **/
   );


/** interface method to detect the structure including presolving
 * @returns SCIP return code */
extern
SCIP_RETCODE DECdetectStructure(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_RESULT*          result             /**< Result pointer to indicate whether some structure was found */
   );


/** writes all selected decompositions */
extern
SCIP_RETCODE DECwriteSelectedDecomps(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 directory,          /**< directory for decompositions */
   char*                 extension           /**< extension for decompositions */
   );


/** write out all known decompositions
 * @returns SCIP return code  */
extern
SCIP_RETCODE DECwriteAllDecomps(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 directory,          /**< directory for decompositions */
   char*                 extension,          /**< the file extension for the export */
   SCIP_Bool             original,           /**< should decomps for original problem be written */
   SCIP_Bool             presolved           /**< should decomps for preoslved problem be written */
   );


/**
 * method to eliminate duplicate constraint names and name unnamed constraints
 * @return SCIP return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompRepairConsNames(
   SCIP*                scip  /**< SCIP data structure */
   );


/** gets an array of all seeeds that are currently considered relevant
 * @params seeedswr  output of the relevant seeeds (don't forget to free the individual wrappers after use)
 * @params nseeeds   amount of seeeds that are put in the array
 * @retruns SCIP return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompGetAllRelevantSeeeds(
   SCIP* scip,                /**< SCIP data structure */
   SEEED_WRAPPER** seeedswr,  /**< seeed wrapper array for output */
   int* nseeeds               /**< number of seeeds in output */
   );

/** write the seeed in conshdlrdata->seeedtowrite into the file as dec
 *
 * @returns SCIP return code */
extern
SCIP_RETCODE SCIPconshdlrDecompWriteDec(
   SCIP*     scip,         /**< SCIP data structure */
   FILE*     file,         /**< file for output */
   SCIP_Bool transformed,  /**< is the problem transformed yet */
   SCIP_RESULT* result     /**< result of writing dec */
   );


/** Runs the bender detector to create a block matrix and outputs its visualization as .png file
 * @returns SCIP return code*/
extern
SCIP_RETCODE SCIPconshdlrDecompWriteMatrix(
   SCIP*                 scip,               /**< scip data structure */
   const char*           filename,           /**< filename the output should be written to (including directory) */
   const char*           workfolder,         /**< directory in which should be worked */
   SCIP_Bool             originalmatrix      /**< should the original (or transformed) matrix be written */
);


/** Gets the best known decomposition
 * @note caller has to free returned DEC_DECOMP
 * @returns the decomposition if available and NULL otherwise */
extern
DEC_DECOMP* DECgetBestDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** Gets the currently considered best seeed
 * @returns the Seeed of the best Seeed if available and seeedwrapper->seeed = NULL otherwise */
extern
SCIP_RETCODE DECgetSeeedToWrite(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             transformed,        /**< is the problem transformed yet */
   SEEED_WRAPPER*        seeedwrapper        /**< seeed wrapper to output */
   );

/** writes out a list of all detectors
 * @returns nothing */
extern
void DECprintListOfDetectors(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets all currently finished decomps
 * @note: the array is allocated and needs to be freed after use!
 * @returns an array containg all finished decompositions
 * */
extern
DEC_DECOMP** SCIPconshdlrDecompGetFinishedDecomps(
   SCIP*     scip /**< SCIP data structure */
   );

/** Gets the number of all finished Seeeds
 * @returns number of finished Seeeds */
extern
int SCIPconshdlrDecompGetNFinishedDecomps(
   SCIP*       scip  /**< SCIP data structure */
   );

/** Gets the number of all seeeds
 * @returns number of Seeeds */
extern
int SCIPconshdlrDecompGetNSeeeds(
   SCIP*       scip  /**< SCIP data structure */
   );

/** Gets the number of all detectors
 * @returns number of detectors */
extern
int SCIPconshdlrDecompGetNDetectors(
   SCIP* scip  /**< SCIP data structure */
   );

/** Gets whether the detection already took place
 * @returns true if detection took place, false otherwise */
extern
SCIP_Bool GCGdetectionTookPlace(
   SCIP*  scip /**< SCIP data structure */
     );

/** Gets an array of all detectors
 *
 * @returns array of detectors */
extern
DEC_DETECTOR** SCIPconshdlrDecompGetDetectors(
   SCIP* scip  /**< SCIP data structure */
   );


/** Gets whether the detection has been performed
 * @returns whether the detection has been performed */
extern
SCIP_Bool DEChasDetectionRun(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** Gets the character of the detector
 * @returns detector character */
extern
char DECdetectorGetChar(
   DEC_DETECTOR*         detector            /**< pointer to detector */
);

/** display statistics about detectors
 * @returns SCIP return code */
extern
SCIP_RETCODE GCGprintDetectorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
   );

/** sets detector parameters values to
 *
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all detector parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for detection is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the detectors produce more decompositions
 *  - SCIP_PARAMSETTING_OFF which turns off all detection
 *
 *   @returns SCIP return code
 */
extern
SCIP_RETCODE GCGsetDetection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );


/** Gets wrapped Seeed with given id
 *
 * @returns SCIP return code */
extern
SCIP_RETCODE GCGgetSeeedFromID(
   SCIP*          scip,       /**< SCIP data structure */
   int*           seeedid,    /**< id of Seeed */
   SEEED_WRAPPER* seeedwr     /**< wrapper for output Seeed */
   );


/* public methods for internal management of seeeds in explore menu and related functions */

/**
 * @brief initilizes the candidates data structures with selected seeeds (or all if there are no selected seeeds) and sort them according to the current scoretype
 * @param scip SCIP data structure
 * @param updatelist whether or not the seeed list should be updated
 * @returns SCIP return code
 */
 extern
SCIP_RETCODE SCIPconshdlrDecompChooseCandidatesFromSelected(
   SCIP* scip,
   SCIP_Bool updatelist
   );

/**
 * @brief method to update the list of incomplete decompositions
 *
 * this list changes due to new decompositions, modified, decompositions or changes of the score
 * @param scip SCIP data structure
 * @returns SCIP return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompUpdateSeeedlist(
   SCIP*                 scip
   );

/**
 * Gets the currently selected scoretype
 * @returns the currently selected scoretype
 */
SCORETYPE SCIPconshdlrDecompGetScoretype(
   SCIP*          scip  /**< SCIP data structure */
   );

/** @brief Gets the first id to visualize in explore menu
 *  @returns id to start with */
int GCGgetSelectFirstIdToVisu(
   SCIP*          scip  /**< SCIP data structure */
   );


/** @brief sets the first id to visualize in explore menu */
void GCGsetSelectFirstIdToVisu(
   SCIP*          scip,  /**< SCIP data structure */
   int            id     /**< id to start at */
   );

/** @brief Gets number of decompositions to be displayed at once in explore menu
 *  @returns id to start with */
int GCGgetSelectVisuLength(
   SCIP*          scip  /**< SCIP data structure */
   );


/** @brief sets number of decompositions to be displayed at once in explore menu */
void GCGsetSelectVisuLength(
   SCIP*          scip,    /**< SCIP data structure */
   int            length   /**< number of rows */
   );

/** @brief Gets a list of ids of the current decomps to be shown in explore menu
 *
 * (to visualize, write, consider for family tree, consider for solving etc. )
 *  @returns list of seeeds */
int** GCGgetSelectList(
   SCIP*          scip  /**< SCIP data structure */
   );


/** @brief Gets a vector containing the indices of selected decompositions in explore menu
 *  @returns list of seeeds */
std::vector<int>* GCGgetSelectIds(
   SCIP*          scip  /**< SCIP data structure */
   );


/** @brief Sets the vector containing the indices of selected decompositions in explore menu */
void GCGsetSelectIds(
   SCIP*          scip,          /**< SCIP data structure */
   std::vector<int>* list   /**< current list of seeeds */
   );

/** @brief Gets whether there are selected decompositions
 *  @returns true iff there are selected decompositions */
SCIP_Bool GCGgetSelectExists(
   SCIP*          scip  /**< SCIP data structure */
   );


/** @brief Sets whether there are selected decompositions */
void GCGsetSelectExists(
   SCIP*          scip,          /**< SCIP data structure */
   SCIP_Bool      selected       /**< input true iff there are selected decompositions */
   );


/** @brief Gets current user seeed
 *
 * (seeed currently selected for modification in explore menu)
 *  @returns current user seeed */
SeeedPtr GCGgetSelectCurrUserSeeed(
   SCIP*          scip  /**< SCIP data structure */
   );

/** @brief Sets current user seeed
 *
 * (seeed currently selected for modification in explore menu)
 *  @returns current user seeed */
void GCGsetSelectCurrUserSeeed(
   SCIP*          scip,    /**< SCIP data structure */
   SeeedPtr       seeed    /**< current seeed */
   );

/** @brief Gets last user seeed
 *
 * (last seeed selected for modification in explore menu)
 *  @returns last user seeed */
SeeedPtr GCGgetSelectLastUserSeeed(
   SCIP*          scip  /**< SCIP data structure */
   );

/** @brief Sets last user seeed
 *
 * (last seeed selected for modification in explore menu)
 *  @returns last user seeed */
void GCGsetSelectLastUserSeeed(
   SCIP*          scip,    /**< SCIP data structure */
   SeeedPtr       seeed    /**< current seeed */
   );

#ifdef __cplusplus
}
#endif

#endif
