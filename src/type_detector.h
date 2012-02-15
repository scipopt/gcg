/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_detector.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for branching rules in gcg projects
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_DETECTOR_H__
#define __SCIP_TYPE_DETECTOR_H__

#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_result.h"
#include "type_decomp.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DEC_Detector DEC_DETECTOR;
typedef struct DEC_DetectorData DEC_DETECTORDATA;

/**
 * initialize data for a detector
 *
 * scip SCIP data structure
 */
#define DEC_DECL_INITDETECTOR(x) SCIP_RETCODE x (SCIP* scip, DEC_DETECTOR* detector)

/**
 * set decomp structure
 */
#define DEC_DECL_SETSTRUCTDECOMP(x) void x (SCIP* scip, DECDECOMP* decdecomp)

/** free data from a detector */
#define DEC_DECL_EXITDETECTOR(x) SCIP_RETCODE x (SCIP* scip, DEC_DETECTOR* detector)

/** detects the structure of a the problem */
#define DEC_DECL_DETECTSTRUCTURE(x) SCIP_RETCODE x (SCIP* scip, DEC_DETECTORDATA* detectordata, DECDECOMP*** decdecomps, int* ndecdecomps, SCIP_RESULT* result)

/** get detector priority */
#define DEC_DECL_GETPRIORITY(x) int x (SCIP* scip, DEC_DETECTORDATA* detectordata)

/** get detector enable/disable status */
#define DEC_DECL_GETISENABLED(x) int x (SCIP* scip, DEC_DETECTORDATA* detectordata)

#ifdef __cplusplus
}
#endif

#endif
