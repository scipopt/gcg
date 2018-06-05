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


/*!
 * \brief help enum to avoid code duplication for the toolbox methods of the detectors
 */
enum toolboxtype {
   PROPAGATE,
   FINISH,
   POSTPROCESS
};



/** forward declarations */
struct seeedpool_wrapper;
typedef struct seeedpool_wrapper SEEEDPOOL_WRAPPER ;

struct Seeed_Wrapper;
typedef struct Seeed_Wrapper SEEED_WRAPPER;


/*!
 * \brief sort the finished decompositions according to the currently chosen score in the according datastructures for the presolved and original problem
 */
SCIP_RETCODE DECconshdlrDecompSortDecompositionsByScore(
   SCIP*          scip
);


/** creates the constraint handler for decomp and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the decomposition structure **/
extern
DEC_DECOMP** SCIPconshdlrDecompGetDecdecomps(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the decomposition structure **/
extern
int SCIPconshdlrDecompGetNDecdecomps(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of conss that were active while detecting decomp originating from seeed with given id **/
extern
int SCIPconshdlrDecompGetNFormerDetectionConssForID(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   id                  /**< id of the seeed */
   );

extern
const char* SCIPconshdlrDecompGetPdfReader(
   SCIP*                scip
   );

extern
SCIP_RETCODE SCIPconshdlrDecompNotifyNonFinalFreeTransform(
   SCIP*                scip
   );

extern
SCIP_RETCODE SCIPconshdlrDecompNotifyFinishedNonFinalFreeTransform(
   SCIP*                scip
   );


/** returns the data of the provided detector */
extern
DEC_DETECTORDATA* DECdetectorGetData(
   DEC_DETECTOR*         detector            /**< Detector data structure */
   );

/** returns the name of the provided detector */
extern
const char* DECdetectorGetName(
   DEC_DETECTOR*         detector            /**< detector data structure */
   );

/** searches for the detector and returns it or returns NULL if detector is not found*/
extern
DEC_DETECTOR* DECfindDetector(
   SCIP*                 scip,               /**< SCIP data structure  */
   const char*           name                /**< name of the detector */
   );

/** includes the detector */
extern
SCIP_RETCODE DECincludeDetector(
   SCIP*                 scip,               /**< SCIP data structure                                                */
   const char*           name,               /**< name of the detector                                               */
   const char            decchar,            /**< display character of the detector                                  */
   const char*           description,        /**< description of the detector                                        */
   int                   freqCallRound,      /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
   int                   maxCallRound,       /** last round the detector gets called                              */
   int                   minCallRound,       /** first round the detector gets called (offset in detection loop) */
   int                   freqCallRoundOriginal,  /** frequency the detector gets called in detection loop while detecting of the original problem */
   int                   maxCallRoundOriginal,   /** last round the detector gets called while detecting of the original problem */
   int                   minCallRoundOriginal,   /** first round the detector gets called (offset in detection loop) while detecting of the original problem */
   int                   priority,           /**< priority of the detector                                           */
   SCIP_Bool             enabled,            /**< whether the detector should be enabled by default                  */
   SCIP_Bool             enabledOriginal,        /**< whether the detector should be enabled by default for detecting the original problem */
   SCIP_Bool             enabledFinishing,   /**< whether the finishing should be enabled */
   SCIP_Bool             enabledPostprocessing,   /**< whether the postprocessing should be enabled */
   SCIP_Bool             skip,               /**< whether the detector should be skipped if others found structure   */
   SCIP_Bool             usefulRecall,       /** is it useful to call this detector on a descendant of the propagated seeed */
   SCIP_Bool             legacymode,         /**< whether (old) DETECTSTRUCTURE method should also be used for detection */
   DEC_DETECTORDATA      *detectordata,      /**< the associated detector data (or NULL)                             */
   DEC_DECL_DETECTSTRUCTURE((*detectStructure)), /**< the method that will detect the structure (must not be NULL)   */
   DEC_DECL_FREEDETECTOR((*freeDetector)),   /**< destructor of detector (or NULL) */
   DEC_DECL_INITDETECTOR((*initDetector)),   /**< initialization method of detector (or NULL)                        */
   DEC_DECL_EXITDETECTOR((*exitDetector)),    /**< deinitialization method of detector (or NULL)                      */
   DEC_DECL_PROPAGATESEEED((*propagateSeeedDetector)),
   DEC_DECL_PROPAGATEFROMTOOLBOX((*propagateFromToolboxDetector)),
   DEC_DECL_FINISHFROMTOOLBOX((*finishFromToolboxDetector)),
   DEC_DECL_FINISHSEEED((*finishSeeedDetector)),
   DEC_DECL_POSTPROCESSSEEED((*postprocessSeeedDetector)),
   DEC_DECL_SETPARAMAGGRESSIVE((*setParamAggressiveDetector)),
   DEC_DECL_SETPARAMDEFAULT((*setParamDefaultDetector)),
   DEC_DECL_SETPARAMFAST((*setParamFastDetector))
   );

/** returns the remaning time of scip that the decomposition may use */
extern
SCIP_Real DECgetRemainingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );



/** checks if two pricing problems are identical based on information from detection */
extern
SCIP_RETCODE SCIPconshdlrDecompArePricingprobsIdenticalForSeeedid(
   SCIP*                scip,
   int                  seeedid,
   int                  probnr1,
   int                  probnr2,
   SCIP_Bool*           identical
   );

/** for two identical pricing problems a corresponding varmap is created */
extern
SCIP_RETCODE SCIPconshdlrDecompCreateVarmapForSeeedId(
   SCIP*                scip,
   SCIP_HASHMAP**       hashorig2pricingvar, /**< mapping from orig to pricingvar  */
   int                  seeedid,
   int                  probnr1,
   int                  probnr2,
   SCIP*                scip1,
   SCIP*                scip2,
   SCIP_HASHMAP*        varmap
   );


/** sets (and adds) the decomposition structure **/
extern
SCIP_RETCODE SCIPconshdlrDecompAddDecdecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   );


/** creates the seeedpool for the presolved problem **/
SCIP_RETCODE SCIPconshdlrDecompCreateSeeedpool(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** creates the seeedpool for the unpresolved problem **/
SCIP_RETCODE SCIPconshdlrDecompCreateSeeedpoolUnpresolved(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** @TODO: IMPORTANT this method will be deleted if the corresponidng wrapper classes are introduced **/
SEEEDPOOL_WRAPPER* SCIPconshdlrDecompGetSeeedpoolUnpresolvedExtern(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** @TODO: IMPORTANT this method will be deleted if the corresponidng wrapper classes are introduced **/
SEEEDPOOL_WRAPPER* SCIPconshdlrDecompGetSeeedpoolExtern(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** creates a user seeed for the presolved problem **/
SCIP_RETCODE SCIPconshdlrDecompCreateUserSeeed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             presolved,          /**< should the user seeed be created for the presolved problem */
   SCIP_Bool             markedincomplete
   );

SCIP_RETCODE SCIPconshdlrDecompExecSelect(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
   );

SCIP_RETCODE SCIPconshdlrDecompExecToolboxModify(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
   );

SCIP_RETCODE SCIPconshdlrDecompExecToolboxCreate(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
   );

SCIP_RETCODE SCIPconshdlrDecompExecToolbox(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
   );


SCIP_Bool SCIPconshdlrDecompUnpresolvedUserSeeedAdded(
   SCIP*                 scip                /**< SCIP data structure */
   );

SCIP_Bool SCIPconshdlrDecompUnpresolvedSeeedExists(
   SCIP*                 scip                /**< SCIP data structure */
   );


SCIP_RETCODE   SCIPconshdlrDecompPopulateSelected(
   SCIP*       scip
   );

SCIP_RETCODE SCIPconshdlrDecompUpdateSeeedlist(
   SCIP*                 scip
   );

SCIP_Bool SCIPconshdlrDecompHasDecomp(
   SCIP*    scip
   );

/** sets the number of blocks */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetnumberOfBlocks(
   SCIP*                 scip,                /**< SCIP data structure */
   int                   nblocks              /**< number of blocks */
   );

/** returns whether there is an user seeed  */
SCIP_Bool SCIPconshdlrDecompUserSeeedIsActive(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** sets the number of blocks */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsDefaultMaster(
   SCIP*                 scip,                /**< SCIP data structure */
   SCIP_Bool             consdefaulttomaster  /**< are not specified constraints set to master for default */
   );



/** sets a constraint by name to a block in the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsToBlock(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           consname,            /**< name of the constraint */
   int                   blockid              /* block index ( counting from 0) */
   );

/** sets a constraint by name to the master in the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsToMaster(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           consname
   );

/** sets a variable by name to a block in the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToBlock(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           varname,             /**< name of the variable */
   int                   blockid              /**< block index ( counting from 0) */
   );

/** sets a variable by name to the master in the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToMaster(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           varname              /**< name of the variable */
   );

/** sets a variable by name to the linking variables in the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToLinking(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           varname              /**< name of the variable */
   );


SCIP_RETCODE SCIPconshdlrDecompAddBlockNumberCandidate(
   SCIP*                 scip,                /**< SCIP data structure */
   int                   blockNumberCandidate /**< name of the variable */
   );

 int SCIPconshdlrDecompGetNBlockNumberCandidates(
   SCIP*                 scip                /**< SCIP data structure */
    );

 int SCIPconshdlrDecompGetBlockNumberCandidate(
    SCIP*                 scip,                /**< SCIP data structure */
    int                   index
     );

 extern
 SCIP_Real SCIPconshdlrDecompGetCompleteDetectionTime(
    SCIP*                 scip
    );

/** rejects and deletes the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedReject(
   SCIP*                 scip                 /**< SCIP data structure */
   );


/** finalizes and flushes the current user seeed, i.e. consider implicits, calc hashvalue, construct decdecomp if complete etc */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedFlush(
   SCIP*                 scip                 /**< SCIP data structure */
   );


SCIP_RETCODE SCIPconshdlrDecompTranslateAndAddCompleteUnpresolvedSeeeds(
   SCIP*                 scip,                 /**< SCIP data structure */
   SCIP_Bool*            success
   );


SCIP_Bool SCIPconshdlrDecompExistsSelected(
   SCIP* scip
   );


int SCIPconshdlrDecompIncreaseAndGetNCallsCreateDecomp(
  SCIP*                 scip                /**< SCIP data structure **/
   );

int SCIPconshdlrDecompDecreaseAndGetNCallsCreateDecomp(
  SCIP*                 scip                /**< SCIP data structure **/
   );


SCIP_RETCODE SCIPconshdlrDecompChooseCandidatesFromSelected(
   SCIP* scip,
   SCIP_Bool updatelist
   );

/** calls old detectStructure methods of chosen detectors, translates the resulting decompositions
 *  into seeeds and adds these seeeds to (presolved) seeedpool */
SCIP_RETCODE SCIPconshdlrDecompAddLegacymodeDecompositions(
   SCIP* scip,
   SCIP_RESULT result
   );

SCIP_Bool SCIPconshdlrDecompDetectBenders(
   SCIP*                   scip
   );


SCIP_Bool SCIPconshdlrDecompIsBestCandidateUnpresolved(
   SCIP*                   scip
   );


SCIP_Bool SCIPconshdlrDecompCheckConsistency(
   SCIP* scip
   );

/** returns the next seeed id managed by cons_decomp */
int SCIPconshdlrDecompGetNextSeeedID(
   SCIP*   scip
   );


SCORETYPE SCIPconshdlrDecompGetCurrScoretype(
   SCIP* scip
   );

//    char*  SCIPconshdlrDecompGetScoretypeShortName(
//      SCIP*       scip,
//      SCORETYPE   sctype
//      );
//
//   char*  SCIPconshdlrDecompGetScoretypeDescription(
//      SCIP*          scip,
//      SCORETYPE      sctype
//         );



/** interface method to detect the structure */
extern
SCIP_RETCODE DECdetectStructure(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result             /**< Result pointer to indicate whether some structure was found */
   );


/** write out all known decompositions **/
SCIP_RETCODE DECwriteAllDecomps(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 directory,          /**< directory for decompositions */
   char*                 extension,          /**< the file extension for the export */
   SCIP_Bool             original,           /**< should decomps for original problem be written */
   SCIP_Bool             presolved           /**< should decomps for preoslved problem be written */

   );


/** gets an array of all seeeds that are currently considered relevant
 * @params seeedswr  output of the relevant seeeds (don't forget to free the individual wrappers after use)
 * @params nseeeds   amount of seeeds that are put in the array
 */
SCIP_RETCODE SCIPconshdlrDecompGetAllRelevantSeeeds(
   SCIP* scip,                /**< SCIP data structure */
   SEEED_WRAPPER** seeedswr,  /**< seeed wrapper array for output */
   int* nseeeds               /**< number of seeeds in output */
   );


/** write family tree **/
SCIP_RETCODE DECwriteFamilyTree(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< filename the output should be written to (including directory) */
   const char*           workfolder,         /**< directory in which should be worked */
   int                   ndecompositions,    /**< the number of (complete) decompositions in order of a certain measure (atm: max white) */
   SCIP_Bool 			    draft               /**< draft mode will not visualize non-zeros but is faster and takes less memory */
   );

extern
SCIP_RETCODE SCIPconshdlrDecompWriteDec(
   SCIP*     scip,
   FILE*     file,
   SCIP_Bool transformed,
   SCIP_RESULT* result
   );

/** returns the best known decomposition, if available and NULL otherwise */
extern
DEC_DECOMP* DECgetBestDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the Seeed ID of the best Seeed if available and -1 otherwise */
SCIP_RETCODE DECgetSeeedToWrite(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             transformed,
   SEEED_WRAPPER*        seeedwrapper        /**< seeed wrapper to output */
   );

/** writes out a list of all detectors */
void DECprintListOfDetectors(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets all currently finished decomps
 * Note: the array is allocated and needs to be freed after use!
 * */
DEC_DECOMP** SCIPconshdlrDecompGetFinishedDecomps(
   SCIP*     scip
   );

/* returns number of finished Seeeds */
int SCIPconshdlrDecompGetNFinishedDecomps(
   SCIP*       scip
   );

/* returns number of all Seeeds */
int SCIPconshdlrDecompGetNSeeeds(
   SCIP*       scip
   );

int SCIPconshdlrDecompGetNDetectors(
   SCIP* scip
   );

SCIP_Bool GCGdetectionTookPlace(
   SCIP*  scip
     );


DEC_DETECTOR** SCIPconshdlrDecompGetDetectors(
   SCIP* scip
   );


/** returns whether the detection has been performed */
SCIP_Bool DEChasDetectionRun(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the character of the detector */
char DECdetectorGetChar(
   DEC_DETECTOR*         detector            /**< pointer to detector */
);

/** display statistics about detectors */
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
 */
SCIP_RETCODE GCGsetDetection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );


/** returns wrapped Seeed with given id */
SCIP_RETCODE GCGgetSeeedFromID(
   SCIP*          scip,       /**< SCIP data structure */
   int*           seeedid,    /**< id of Seeed */
   SEEED_WRAPPER* seeedwr     /**< wrapper for output Seeed */
   );


/** returns wrapped Seeedpools */
SCIP_RETCODE GCGgetCurrentSeeedpools(
   SCIP*          scip,                   /**< SCIP data structure */
   SEEED_WRAPPER* seeedpoolwr,            /**< wrapper for presolved output Seeedpool (or NULL) */
   SEEED_WRAPPER* seeedpoolunpresolvedwr  /**< wrapper for unpresolved output Seeedpool (or NULL) */
   );


/** gets the ids of all selected seeeds */
SCIP_RETCODE SCIPconshdlrDecompGetSelectedSeeeds(
   SCIP* scip,       /**< SCIP data structure */
   int** output,     /**< array to put ids into */
   int* outputsize   /**< size of output */
   );

SCIP_RETCODE GCGprintMiplibBaseInformation(
   SCIP*                scip,
   FILE*                file
   );

SCIP_RETCODE GCGprintMiplibBaseInformationHeader(
   SCIP*                scip,
   FILE*                file
   );


SCIP_RETCODE GCGprintMiplibConnectedInformation(
   SCIP*                scip,
   FILE*                file
   );

SCIP_RETCODE GCGprintMiplibDecompInformation(
   SCIP*                scip,
   FILE*                file
   );

SCIP_RETCODE GCGprintOptionalOutput(
   SCIP*                scip,
   SCIP_DIALOGHDLR*     dialoghdlr         /**< dialog handler */
   );


#ifdef __cplusplus
}
#endif

#endif
