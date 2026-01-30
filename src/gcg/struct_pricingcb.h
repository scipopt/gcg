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

/**@file   struct_pricingcb.h
 * @ingroup INTERNALAPI
 * @brief  data structure for pricing callbacks
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_PRICINGCB_H__
#define __SCIP_STRUCT_PRICINGCB_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "gcg/type_pricingcb.h"

#ifdef __cplusplus
extern "C" {
#endif

/** pricing callbacks data */
struct GCG_Pricingcb
{
   SCIP_Longint          nprepricingcalls;   /**< number of times, the pre-pricing method was called */
   SCIP_Longint          npostpricingcalls;  /**< number of times, the post-pricing method was called */
   char*                 name;               /**< name of the pricing callback */
   char*                 desc;               /**< description of the pricing callback */
   GCG_DECL_PRICINGCBCOPY((*pricingcbcopy));/**< copy method of the pricing callback or NULL if you don't want to copy your plugin into sub-SCIPs */
   GCG_DECL_PRICINGCBFREE((*pricingcbfree));/**< destructor of the pricing callback */
   GCG_DECL_PRICINGCBINIT((*pricingcbinit));/**< initialize the pricing callback */
   GCG_DECL_PRICINGCBEXIT((*pricingcbexit));/**< deinitialize the pricing callback */
   GCG_DECL_PRICINGCBINITSOL((*pricingcbinitsol));/**< solving process initialization method of the pricing callback */
   GCG_DECL_PRICINGCBEXITSOL((*pricingcbexitsol));/**< solving process deinitialization method of the pricing callback */
   GCG_DECL_PRICINGCBPREPRICING((*pricingcbprepricing));/**< the pre-pricing method of the pricing callback */
   GCG_DECL_PRICINGCBPOSTPRICING((*pricingcbpostpricing));/**< the post-pricing method of the pricing callback */
   GCG_PRICINGCBDATA*   pricingcbdata;      /**< pricing callbacks local data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up the pricing callback plugin */
   SCIP_CLOCK*           pricingcbclock;     /**< the execution time of the pricing callback plugin */
   int                   priority;           /**< priority of the pricing callbacks */

   /* additional pricing callback parameters */
   SCIP_Bool             enabled;            /**< is this pricing callback enabled? */
   SCIP_Bool             exclusive;          /**< is this pricing callback executed exclusively? (must be highest priority) */
};

#ifdef __cplusplus
}
#endif

#endif
