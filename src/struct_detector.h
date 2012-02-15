/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_detector.h
 * @brief  data structures for detectors
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_DETECTOR_H__
#define __SCIP_STRUCT_DETECTOR_H__

#include "type_detector.h"

#ifdef __cplusplus
extern "C" {
#endif

/** detector data structure
 *
 * @todo: add priority and enabled flag (bug #7)
 */
struct DEC_Detector {

   DEC_DECL_INITDETECTOR((*initDetection));     /**< initialization method of detector */
   DEC_DECL_DETECTSTRUCTURE((*detectStructure));/**< structure detection method of detector */
   DEC_DECL_EXITDETECTOR((*exitDetection));     /**< deinitialization method of detector */
   DEC_DECL_GETPRIORITY((*getPriority));        /**< getter method for the priority */
   DEC_DECL_GETISENABLED((*getIsEnabled));      /**< getter method to check whether the detector is enabled */
   DEC_DETECTORDATA* decdata;                   /**< custom data structure of the detectors */
   const char *name;                            /**< name of the detector */
   char decchar;                                /**< display character of detector */
};


#ifdef __cplusplus
}
#endif

#endif
