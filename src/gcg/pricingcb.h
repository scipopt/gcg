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

/**@file   pricingcb.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for the pricing callback
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRICINGCB_H__
#define __SCIP_PRICINGCB_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "gcg/type_pricingcb.h"
#include "gcg/pub_pricingcb.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates a pricing callback */
SCIP_RETCODE GCGpricingcbCreate(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB**       pricingcb,          /**< pointer to the pricing callback data structure */
   const char*           name,               /**< name of the pricing callback */
   const char*           desc,               /**< description of the pricing callback */
   int                   priority,           /**< priority of the pricing callback */
   GCG_DECL_PRICINGCBFREE((*pricingcbfree)), /**< destructor of the pricing callback */
   GCG_DECL_PRICINGCBINIT((*pricingcbinit)), /**< initialize pricing callback */
   GCG_DECL_PRICINGCBEXIT((*pricingcbexit)), /**< deinitialize pricing callback */
   GCG_DECL_PRICINGCBINITSOL((*pricingcbinitsol)),/**< solving process initialization method of the pricing callback */
   GCG_DECL_PRICINGCBEXITSOL((*pricingcbexitsol)),/**< solving process deinitialization method of the pricing callback */
   GCG_DECL_PRICINGCBPREPRICING((*pricingcbprepricing)),/**< pre-pricing method of the pricing callback */
   GCG_DECL_PRICINGCBPOSTPRICING((*pricingcbpostpricing)),/**< post-pricing method of the pricing callback */
   GCG_PRICINGCBDATA*  pricingcbdata         /**< pricing callback data */
   );

/** calls destructor and frees memory of the pricing callback */
SCIP_RETCODE GCGpricingcbFree(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB**       pricingcb           /**< pointer to the pricing callback data structure */
   );

/** initializes the pricing callback */
SCIP_RETCODE GCGpricingcbInit(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb           /**< the pricing callback */
   );

/** calls exit method of the pricing callback */
SCIP_RETCODE GCGpricingcbExit(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );

/** informs the pricing callback that the branch and bound process is being started */
SCIP_RETCODE GCGpricingcbInitsol(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );

/** informs the pricing callback that the branch and bound process data is being freed */
SCIP_RETCODE GCGpricingcbExitsol(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb           /**< pricing callback */
   );

/** calls pre-pricing method of the pricing callback */
SCIP_RETCODE GCGpricingcbPrepricing(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   SCIP_PRICER*          pricer,             /**< the pricer */
   GCG_PRICETYPE         type,               /**< the type of pricing, either redcost or farkas */
   SCIP_Bool*            abort,              /**< should the pricing be aborted? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls post-pricing method of the pricing callback */
SCIP_RETCODE GCGpricingcbPostpricing(
   GCG*                  gcg,                /**< GCG data structure */
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   SCIP_PRICER*          pricer,             /**< the pricer */
   GCG_PRICETYPE         type,               /**< the type of pricing, either redcost or farkas */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of the pricing callback */
void GCGpricingcbSetPriority(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   int                   priority            /**< new priority of the pricing callback */
   );

/** sets destructor callback of the pricing callback */
void GCGpricingcbSetFree(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   GCG_DECL_PRICINGCBFREE((*pricingcbfree))  /**< destructor of the pricing callback */
   );

/** sets initialization callback of the pricing callback */
void GCGpricingcbSetInit(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   GCG_DECL_PRICINGCBINIT((*pricingcbinit))  /**< initialize the pricing callback */
   );

/** sets deinitialization callback of the pricing callback */
void GCGpricingcbSetExit(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   GCG_DECL_PRICINGCBEXIT((*pricingcbexit))  /**< deinitialize the pricing callback */
   );

/** sets solving process initialization callback of the pricing callback */
void GCGpricingcbSetInitsol(
   GCG_PRICINGCB*        pricingcb,          /**< pricing callback */
   GCG_DECL_PRICINGCBINITSOL((*pricingcbinitsol))/**< solving process initialization callback of the pricing callback */
   );

/** sets solving process deinitialization callback of the pricing callback */
void GCGpricingcbSetExitsol(
   GCG_PRICINGCB*        pricingcb,          /**< the pricing callback */
   GCG_DECL_PRICINGCBEXITSOL((*pricingcbexitsol))/**< solving process deinitialization callback of the pricing callback */
   );

#ifdef __cplusplus
}
#endif

#endif
