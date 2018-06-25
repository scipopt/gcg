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


/**
 * \brief sort the finished decompositions according to the currently chosen score in the according datastructures for the presolved and original problem
 * @param scip data stucture
 * @return scip return code
 */
/*!
 *
 */
SCIP_RETCODE DECconshdlrDecompSortDecompositionsByScore(
   SCIP*          scip
);


/**
 * @brief creates the constraint handler for decomp and includes it in SCIP
 * @param scip scip date structure
 * @return scip return code
 */
extern
SCIP_RETCODE SCIPincludeConshdlrDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**
 * returns the decomposition structures, caller is responsible for freeing memory
 * @param scip scip data structure
 * @return scip return code
 */
extern
DEC_DECOMP** SCIPconshdlrDecompGetDecdecomps(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**
 * returns the number of found decomposition structures
 * @param scipscip dat structures
 * @return scip return code
 */
extern
int SCIPconshdlrDecompGetNDecdecomps(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**
 * @brief returns the number of constraints that were active while detecting the decomposition originating from the seeed with the given id, this method is used to decide if the problem has changed since detection, if so the aggregation information needs to be recalculated
 * @param scip scip data structure
 * @param id id of the seeed the information is asked for
 * @return number of constraints that were active while detecting the decomposition originating from the seeed with the given id, this method is used to decide if the problem has changed since detection, if so the aggregation information needs to be recalculated
 */
extern
int SCIPconshdlrDecompGetNFormerDetectionConssForID(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   id                  /**< id of the seeed */
   );

/**
 * @brief returns string name of the chosen pdf reader
 *
 * \ref see parameter visual/pdfreader
 * @param scip scip data structure
 * @return returns string name of the chosen pdf reader
 */
extern
const char* SCIPconshdlrDecompGetPdfReader(
   SCIP*                scip
   );


/**
 * @brief used before calling SCIPfreeTransform(),, if called to revoke presolving (e.g. if unpresolved decomposition is used, and transformation is not successful), this seems mandatory to decide during consExitDecomp if the original detection information should be freed
 * @param scip data structure
 * @return scip return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompNotifyNonFinalFreeTransform(
   SCIP*                scip
   );


/**
 *
 * @brief used after calling SCIPfreeTransform() if called to revoke presolving (e.g. if unpresolved decomposition is used, and transformation is not successful), this seems mandatory to decide during consExitDecomp if the original detection information should be freed
 * @param scip data structure
 * @return scip return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompNotifyFinishedNonFinalFreeTransform(
   SCIP*                scip
   );


/**
 * @brief returns the data of the provided detector
 * @param detector Detector data structure
 * @return data of the provided detector
 */
extern
DEC_DETECTORDATA* DECdetectorGetData(
   DEC_DETECTOR*         detector            /**< Detector data structure */
   );

/**
 * @brief returns the name of the provided detector
 * @param detector detector data structure
 * @return name of the given detector
 */
extern
const char* DECdetectorGetName(
   DEC_DETECTOR*         detector            /**< detector data structure */
   );

/**
 * @brief searches for the detector with the given name and returns it or NULL if detector is not found
 * @param scip scip data strcuture
 * @param name the name of the searched detector
 * @return returns detector pointer or NULL if detector with given name is not found
 */
extern
DEC_DETECTOR* DECfindDetector(
   SCIP*                 scip,               /**< SCIP data structure  */
   const char*           name                /**< name of the detector */
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
 * @param DEC_DECL_DETECTSTRUCTURE((*detectStructure))   the method that will detect the structure (may be NULL), only used in legecy detection mode
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
 * @return scip return code
 */
/**  */
extern
SCIP_RETCODE DECincludeDetector(
   SCIP*                 scip,                             /**< SCIP data structure                                                */
   const char*           name,                             /**< name of the detector                                               */
   const char            decchar,                          /**< display character of the detector                                  */
   const char*           description,                      /**< description of the detector                                        */
   int                   freqCallRound,                    /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
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
   DEC_DECL_PROPAGATESEEED((*propagateSeeedDetector)),
   DEC_DECL_PROPAGATEFROMTOOLBOX((*propagateFromToolboxDetector)),
   DEC_DECL_FINISHFROMTOOLBOX((*finishFromToolboxDetector)),
   DEC_DECL_FINISHSEEED((*finishSeeedDetector)),
   DEC_DECL_POSTPROCESSSEEED((*postprocessSeeedDetector)),
   DEC_DECL_SETPARAMAGGRESSIVE((*setParamAggressiveDetector)),
   DEC_DECL_SETPARAMDEFAULT((*setParamDefaultDetector)),
   DEC_DECL_SETPARAMFAST((*setParamFastDetector))
   );


/**
 * @brief returns the remaining time of scip that the decomposition may use
 * @param scip scip data structure
 * @return remaining  time that the decompositon may use
 */
extern
SCIP_Real DECgetRemainingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );



/**
 * checks if two pricing problems are identical based on information from detection
 * @param scip scip datta structure
 * @param seeedid if of the partial decompostion for which the pricing problems are checked for identity
 * @param probnr1 index of first block to check
 * @param probnr2 index of second block to check
 * @param identical bool pointer to score the result of the check
 * @return scip return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompArePricingprobsIdenticalForSeeedid(
   SCIP*                scip,             /**< scip scip datta structure */
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
 * @return scip return code
 */
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


/**
 * @brief sets (and adds) the decomposition structure
 * @param scip scip datat structure
 * @param decdecomp decomposition data structure
 * @return scip return code
 */
extern
SCIP_RETCODE SCIPconshdlrDecompAddDecdecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   );

/**
 * @brief creates the seeedpool for the presolved problem
 * @param scip scip data structure
 * @return scip return code
 */
SCIP_RETCODE SCIPconshdlrDecompCreateSeeedpool(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**
 * @brief creates the seeedpool for the unpresolved problem
 * @param scip scip data structure
 * @return scip return code
 */
SCIP_RETCODE SCIPconshdlrDecompCreateSeeedpoolUnpresolved(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**
 * @brief help method to access seeedpool for unpresolved problem
 * @TODO: consider deleting this method will be deleted if the corresponidng wrapper classes are introduced
 * @param scip scip data structure
 * @return pointer to seeedpool wrapper data structure
 */
SEEEDPOOL_WRAPPER* SCIPconshdlrDecompGetSeeedpoolUnpresolvedExtern(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**
 * @brief help method to access seeedpool for transformed problem
 * @TODO: consider deleting this method will be deleted if the corresponidng wrapper classes are introduced
 * @param scip scip dat structure
 * @return pointer to seeedpool wrapper data structure
 */
SEEEDPOOL_WRAPPER* SCIPconshdlrDecompGetSeeedpoolExtern(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**
 * @brief creates a user seeed for the problem
 * @param scip scip data structure
 * @param presolved should the user seeed be created for the presolved (transformed) or problem
 * @param markedincomplete boolean to notify that the created user decomposition is partial and should not be completed by assigning open constraints to the master
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompCreateUserSeeed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             presolved,          /**< should the user seeed be created for the presolved problem */
   SCIP_Bool             markedincomplete    /**< boolean to notify that the created user decomposition is partial and should not be completed by assigning open constraints to the master */
   );

/**
 * @brief method too handle user input for "explore" command
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompExecSelect(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
   );

/**
 * @brief method to handle and moderate user input for modifying decompositions
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompExecToolboxModify(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
   );

/**
 * @brief method to handle and moderate user input for creating new decompositions by the user
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @return SCIP return data structure
 */
SCIP_RETCODE SCIPconshdlrDecompExecToolboxCreate(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
   );

/**
 * @brief method to handle and moderate user input for creating new decompositions and modifying existing decompositions by the user
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompExecToolbox(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
   );


/**
 * @brief returns whether or not an unpresolved (untransformed) decompositions was given by the user
 * @param scip SCIP data structure
 * @return true iff an unpresolved (untransformed) decompositions was given by the user
 */
SCIP_Bool SCIPconshdlrDecompUnpresolvedUserSeeedAdded(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**
 * @brief returns whether or not an unpresolved (untransformed) decompositions exists in the data structures
 * @param scip SCIP data structure
 * @return SCIP return code
 */
SCIP_Bool SCIPconshdlrDecompUnpresolvedSeeedExists(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**
 * @brief populate datastructures with incomplete decompositions (e.g. that were given by the user) to complete them during detection loop
 * @param scip SCIP data structure
 * @return SCIP return code
 */
SCIP_RETCODE   SCIPconshdlrDecompPopulateSelected(
   SCIP*       scip
   );

/**
 * @brief method to update the list of incomplete decompositions in "explore" submenue ( this list changes due to new decompositions,  modified, decompositions or changes of the score
 * @param scip SCIP data structure
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUpdateSeeedlist(
   SCIP*                 scip
   );

/**
 * @brief returns whether or not there exists at least one (com[plete or incomplete) decomposition
 * @param scip SCIP data structure
 * @return TRUE if there exists at least one (com[plete or incomplete) decomposition
 */
SCIP_Bool SCIPconshdlrDecompHasDecomp(
   SCIP*    scip
   );

/**
 * @brief set the number of blocks in the current user seeed (which is used for user input (read or modify) )
 * @param scip SCIP data structure
 * @param nblocks number of blocks
 * @return SCIP return code
 */
/** sets the number of blocks */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetnumberOfBlocks(
   SCIP*                 scip,                /**< SCIP data structure */
   int                   nblocks              /**< number of blocks */
   );


/**
 * @brief returns whether there is an user seeed that is currently worked on
 * @param scip SCIP datastructure
 * @return TRUE there is an user seeed that is currently worked on
 */
SCIP_Bool SCIPconshdlrDecompUserSeeedIsActive(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**
 * @brief set the user given information of the current user seeed according consdefaultmaster (if TRUE all open constraints are set to master )
 * @param scip SCIP data structure
 * @param consdefaulttomaster if TRUE all open constraints are set to master,
 * @return SCIP return code
 */
/** sets the number of blocks */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsDefaultMaster(
   SCIP*                 scip,                /**< SCIP data structure */
   SCIP_Bool             consdefaulttomaster  /**< are not specified constraints set to master for default */
   );


/**
 * @brief sets a constraint by name to a block in the current user seeed
 * @param scip SCIP data structure
 * @param consname name of the constraint that should be set to a block
 * @param blockid index of the block the constraint should be assigned to
 * @return SCIP return code
 */
/** sets a constraint by name to a block in the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsToBlock(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           consname,            /**< name of the constraint */
   int                   blockid              /**< block index ( counting from 0) */
   );

/**
 * @brief sets a constraint by name to master in the current user seeed
 * @param scip SCIP data structure
 * @param consname of the constraint that should be set to master
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsToMaster(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           consname
   );

/**
 * @brief sets a variable by name to a block in the current user seeed
 * @param scip SCIP dat structure
 * @param varname name of the variable
 * @param blockid block index ( counting from 0)
 * @return SCIP return code
 */
/** sets a variable by name to a block in the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToBlock(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           varname,             /**< name of the variable */
   int                   blockid              /**< block index ( counting from 0) */
   );


/**
 * @brief sets a variable by name to the master in the current user seeed
 * @param scip SCIP data structure
 * @param varname name of the variable
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToMaster(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           varname              /**< name of the variable */
   );


/**
 * @brief sets a variable by name to the linking variables in the current user seeed
 * @param scip SCIP data structure
 * @param varname name of the variable
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToLinking(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           varname              /**< name of the variable */
   );

/**
 * @brief add block number user candidate (user candidates are prioritized over found ones)
 * @param scip SCIP data structure
 * @param blockNumberCandidate given block number candidate
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompAddBlockNumberCandidate(
   SCIP*                 scip,                /**< SCIP data structure */
   int                   blockNumberCandidate /**< given block number candidate */
   );


/**
 * @brief returns the number of block candidates given by the user
 * @param scip SCIP data structures
 * @return number of block candidates given by the user
 */
 int SCIPconshdlrDecompGetNBlockNumberCandidates(
   SCIP*                 scip                /**< SCIP data structure */
    );

 /**
  * @brief returns block number user candidate with given index
  * @param scip SCIP data structure
  * @param index index of block number user candidate that should be returned
  * @return block number user candidate with given index
  */
 int SCIPconshdlrDecompGetBlockNumberCandidate(
    SCIP*                 scip,                /**< SCIP data structure */
    int                   index
     );


 /**
  * @brief returns the total detection time
  * @param scip SCIP data structure
  * @return total detection time
  */
 extern
 SCIP_Real SCIPconshdlrDecompGetCompleteDetectionTime(
    SCIP*                 scip
    );

 /**
  * @brief rejects and deletes the current user seeed
  * @param scip SCIP data structure
  * @return SCIP return code
  */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedReject(
   SCIP*                 scip                 /**< SCIP data structure */
   );


/**
 * @brief finalizes and flushes the current user seeed, i.e. consider implicits, calc hashvalue, construct decdecomp if complete etc
 * @param scip SCIP data structure
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedFlush(
   SCIP*                 scip                 /**< SCIP data structure */
   );


/**
 * @brief translates unpresolved seeed to a complete presolved one
 * @param scip SCIP data structure
 * @param success  at least one unpresolved seeed could not be translated in a complete presolved one
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompTranslateAndAddCompleteUnpresolvedSeeeds(
   SCIP*                 scip,                 /**< SCIP data structure */
   SCIP_Bool*            success
   );


/**
 * @brief returns if there is a decomposition that is currently selected by the user (done in explore menue)
 * @param scip SCIP data structure
 * @return TRUE if there is a decomposition that is currently selected by the user (done in explore menue)
 */
SCIP_Bool SCIPconshdlrDecompExistsSelected(
   SCIP* scip
   );


/**
 * @brief counts up the counter for created decompositions and returns it
 * @param scip SCIP data structure
 * @return number of created decompositions that was recently increased
 */
int SCIPconshdlrDecompIncreaseAndGetNCallsCreateDecomp(
  SCIP*                 scip                /**< SCIP data structure **/
   );

/**
 * @brief decreases the counter for created decompositions and returns it
 * @param scip SCIP data structure
 * @return number of created decompositions that was recently decreased
 */
int SCIPconshdlrDecompDecreaseAndGetNCallsCreateDecomp(
  SCIP*                 scip                /**< SCIP data structure **/
   );


/**
 * @brief initilizes the candidates data structures with selected seeeds (or all if there are no selected seeeds) and sort them according to the current scoretype
 * @param scip SCIP data structure
 * @param updatelist whether or not the seeed list should be updated
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompChooseCandidatesFromSelected(
   SCIP* scip,
   SCIP_Bool updatelist
   );


/**
 * @brief calls old detectStructure methods of chosen detectors, translates the resulting decompositions
 *  into seeeds and adds these seeeds to (presolved) seeedpool
 * @param scip SCIP data structure
 * @param result was the legacy call successful
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompAddLegacymodeDecompositions(
   SCIP* scip,
   SCIP_RESULT result
   );


/**
 * @brief returns whether
 * @param scip
 * @return
 */
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
