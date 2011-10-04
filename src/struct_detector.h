/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_branchgcg.h 198 2011-01-06 16:58:56Z ggamrath $"

/**@file   struct_detector.h
 * @brief  datastructures for branching rules
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_DETECTOR_H__
#define __SCIP_STRUCT_DETECTOR_H__

#include "type_detector.h"

#ifdef __cplusplus
extern "C" {
#endif

struct DEC_Detector {

   DEC_DECL_INITDETECTOR((*initDetection));
   DEC_DECL_SETSTRUCTDECOMP((*setStructDecomp));
   DEC_DECL_DETECTSTRUCTURE((*detectStructure));
   DEC_DECL_EXITDETECTOR((*exitDetection));
   DEC_DECL_GETPRIORITY((*getPriority));
   DEC_DETECTORDATA* decdata;
   int i;
   const char *name;
};


#ifdef __cplusplus
}
#endif

#endif
