/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_decomp.h
 * @brief  constraint handler for decomp constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_DECOMP_H__
#define __SCIP_CONS_DECOMP_H__

#include "scip/scip.h"
#include "type_detector.h"


#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for decomp constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrDecomp(
   SCIP* scip                 /**< SCIP data structure */
   );

/** creates and captures a decomp constraint */
extern
SCIP_RETCODE SCIPcreateConsDecomp(
   SCIP*       scip,          /**< SCIP data structure */
   SCIP_CONS** cons,          /**< pointer to hold the created constraint */
   const char* name           /**< name of constraint */
   );

/** returns the decomposition structure **/
extern
DECDECOMP** SCIPconshdlrDecompGetDecdecomps(
      SCIP* scip              /**< SCIP data structure */
   );

/** returns the decomposition structure **/
extern
int SCIPconshdlrDecompGetNDecdecomps(
      SCIP* scip              /**< SCIP data structure */
   );

/** returns the data of the provided detector */
extern
DEC_DETECTORDATA* DECdetectorGetData(
   DEC_DETECTOR* detector     /**< Detector data structure */
   );

/** returns the name of the provided detector */
extern
const char* DECdetectorGetName(
   DEC_DETECTOR*  detector    /**< detector data structure */
   );

/** searches for the detector and returns it or returns NULL if detector is not found*/
extern
DEC_DETECTOR* DECfindDetector(
   SCIP*       scip,          /**< SCIP data structure */
   const char* name           /**< name of the detector */
   );

/** includes the detector */
extern
SCIP_RETCODE DECincludeDetector(
   SCIP* scip,                                     /**< SCIP data structure */
   const char* name,                               /**< name of the detector */
   const char decchar,                             /**< display character of the detector */
   const char* description,                        /**< description of the detector */
   int priority,                                   /**< priority of the detector */
   SCIP_Bool enabled,                              /**< whether the detector should be enabled by default */
   DEC_DETECTORDATA *detectordata,                 /**< the associated detector data (or NULL) */
   DEC_DECL_DETECTSTRUCTURE((*detectStructure)),   /**< the method that will detect the structure (must not be NULL)*/
   DEC_DECL_INITDETECTOR((*initDetector)),         /**< initialization method of detector (or NULL) */
   DEC_DECL_EXITDETECTOR((*exitDetector))          /**< deinitialization method of detector (or NULL) */
   );

/** returns the remaning time of scip that the decomposition may use */
extern
SCIP_Real DECgetRemainingTime(
   SCIP* scip                 /**< SCIP data structure */
   );

/** sets (and adds) the decomposition structure **/
extern
SCIP_RETCODE SCIPconshdlrDecompAddDecdecomp(
   SCIP*      scip,           /**< SCIP data structure */
   DECDECOMP* decdecomp       /**< DECDECOMP data structure */
   );

/** interface method to detect the structure */
extern
SCIP_RETCODE DECdetectStructure(
   SCIP* scip                 /**< SCIP data structure */
   );


/** write out all known decompositions **/
SCIP_RETCODE DECwriteAllDecomps(
   SCIP* scip,                /**< SCIP data structure */
   char* extension            /**< the file extension for the export */
   );

/** returns the best known decomposition, if available and NULL otherwise */
extern
DECDECOMP* DECgetBestDecomp(
   SCIP* scip                 /**< SCIP data structure */
   );

/**< Writes out a list of all detectors */
void DECprintListOfDetectors(
   SCIP* scip                 /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
