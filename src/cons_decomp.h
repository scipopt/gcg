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

typedef struct DEC_Detector DEC_DETECTOR;
typedef struct DEC_DetectorData DEC_DETECTORDATA;


#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for decomp constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a decomp constraint */
extern
SCIP_RETCODE SCIPcreateConsDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name                /**< name of constraint */
   );

/** returns the decomposition structure **/
extern
DECDECOMP* SCIPconshdlrDecompGetDecdecomp(
   SCIP *scip                                /**< SCIP data structure */
   );

/** returns the data of the provided detector */
extern
DEC_DETECTORDATA* DECdetectorGetData(
   DEC_DETECTOR* detector                    /**< Detector data structure */
   );

/** returns the name of the provided detector */
extern
const char* DECdetectorGetName(
   DEC_DETECTOR* detector                    /**< Detector data structure */
   );

/** Searches for the detector and returns it or returns NULL if detector is not found*/
extern
DEC_DETECTOR* DECfindDetector(
   SCIP *scip,                               /**< SCIP data structure */
   const char *name                          /**< Name of the detector to return */
   );

/** includes the detector */
extern
SCIP_RETCODE DECincludeDetector(
   SCIP* scip,                                     /**< SCIP data structure */
   const char *name,                               /**< name of the detector */
   DEC_DETECTORDATA *detectordata,                 /**< the associated detector data (or NULL) */
   DEC_DECL_DETECTSTRUCTURE((*detectStructure)),   /**< the method that will detect the structure (must not be NULL)*/
   DEC_DECL_SETSTRUCTDECOMP((*setStructDecomp)),   /**< interface method to tell detector where to store structure information (must not be NULL) */
   DEC_DECL_INITDETECTOR((*initDetector)),         /**< initialization method of detector (or NULL) */
   DEC_DECL_EXITDETECTOR((*exitDetector)),         /**< deinitialization method of detector (or NULL) */
   DEC_DECL_GETPRIORITY((*getPriority))            /**< interface method to get priority of detector (must not be NULL) */
   );

/** returns the remaning time of scip that the decomposition may use */
extern
SCIP_Real DECgetRemainingTime(
   SCIP* scip                    /**< SCIP data structure */
   );

/** converts the structure to the gcg format by setting the appropriate blocks and master constraints */
extern
SCIP_RETCODE DECOMPconvertStructToGCG(
      SCIP*         scip,     /**< SCIP data structure          */
      DECDECOMP*    decdecomp /**< decdecom data structure      */
   );

#ifdef __cplusplus
}
#endif

#endif
