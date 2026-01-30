/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   struct_colpool.h
 * @ingroup DATASTRUCTURES
 * @brief  data structures for storing cols in a col pool
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_COLPOOL_H__
#define __SCIP_STRUCT_COLPOOL_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_misc.h"
#include "gcg/type_gcgcol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** storage for pooled cols */
struct GCG_Colpool
{
   GCG*                  gcg;                /**< GCG data structure */
   SCIP*                 scip;               /**< SCIP data structure */
   SCIP_Longint          nodenr;             /**< node at which columns in colpool respect branching decisions */
   SCIP_Bool             infarkas;           /**< in Farkas pricing? */
   SCIP_Longint          ncalls;             /**< number of times, the colpool was separated */
   SCIP_Longint          ncolsfound;         /**< total number of cols that were separated from the pool */
   SCIP_CLOCK*           poolclock;          /**< pricing time */
   SCIP_HASHTABLE*       hashtable;          /**< hash table to identify already stored cols */
   GCG_COL**             cols;               /**< stored cols of the pool */
   SCIP_Longint          processedlp;        /**< last LP that has been processed for separating the LP */
   SCIP_Longint          processedlpsol;     /**< last LP that has been processed for separating other solutions */
   int                   colssize;           /**< size of cols array */
   int                   ncols;              /**< number of cols stored in the pool */
   int                   agelimit;           /**< maximum age a col can reach before it is deleted from the pool */
   int                   firstunprocessed;   /**< first col that has not been processed in the last LP */
   int                   firstunprocessedsol;/**< first col that has not been processed in the last LP when separating other solutions */
   int                   maxncols;           /**< maximal number of cols stored in the pool at the same time */
   SCIP_Bool             globalcolpool;      /**< is this the global col pool of SCIP? */
};

#ifdef __cplusplus
}
#endif

#endif
