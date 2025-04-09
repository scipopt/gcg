/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file   pub_pricingcb.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for pricing callback plugins
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_PRICINGCB_H__
#define __SCIP_PUB_PRICINGCB_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "gcg/type_pricingcb.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicBenderscutsMethods
 *
 * @{
 */

/** compares two pricing callback plugins w. r. to their priority */
GCG_EXPORT
SCIP_DECL_SORTPTRCOMP(GCGpricingcbComp);

/** comparison method for sorting pricing callback plugins w.r.t. to their name */
GCG_EXPORT
SCIP_DECL_SORTPTRCOMP(GCGpricingcbCompName);

/** gets user data of the pricing callback plugin */
GCG_EXPORT
GCG_PRICINGCBDATA* GCGpricingcbGetData(
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );

/** sets user data of the pricing callback plugin; user has to free old data in advance! */
GCG_EXPORT
void GCGpricingcbSetData(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   GCG_PRICINGCBDATA*    pricingcbdata       /**< new pricing callback plugin user data */
   );

/** gets name of the pricing callback plugin */
GCG_EXPORT
const char* GCGpricingcbGetName(
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );

/** gets description of the pricing callback plugin */
GCG_EXPORT
const char* GCGpricingcbGetDesc(
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );

/** gets priority of the pricing callback plugin */
GCG_EXPORT
int GCGpricingcbGetPriority(
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );

/** gets the number of times the pre-pricing method of the pricing callback plugin was called */
GCG_EXPORT
SCIP_Longint GCGpricingcbGetNPrepricingCalls(
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );

/** gets the number of times the post-pricing method of the pricing callback plugin was called */
GCG_EXPORT
SCIP_Longint GCGpricingcbGetNPostpricingCalls(
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );

/** gets time in seconds used in this pricing callback plugin for setting up for next stages */
GCG_EXPORT
SCIP_Real GCGpricingcbGetSetupTime(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );

/** gets time in seconds used in this pricing callback plugin */
GCG_EXPORT
SCIP_Real GCGpricingcbGetTime(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );

/** sets the enabled flag of the pricing callback plugin method */
GCG_EXPORT
void GCGpricingcbSetEnabled(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   SCIP_Bool             enabled             /**< flag to indicate whether the pricing callback plugin is enabled */
   );

/** sets the exclusive flag of the pricing callback plugin method */
GCG_EXPORT
void GCGpricingcbSetExclusive(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   SCIP_Bool             exclusive           /**< flag to indicate whether the pricing callback plugin is executed exclusively */
   );

/** returns whether the pricing callback is enabled */
GCG_EXPORT
SCIP_Bool GCGpricingcbIsEnabled(
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );

/** returns whether the methods of this pricing callback should be executed exclusively */
GCG_EXPORT
SCIP_Bool GCGpricingcbIsExclusive(
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );


/** @} */

#ifdef __cplusplus
}
#endif

#endif
