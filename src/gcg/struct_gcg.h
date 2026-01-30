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

/**@file   struct_gcg.h
 * @ingroup DATASTRUCTURES
 * @brief  gcg data structure
 * @author Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_GCG_H_
#define GCG_STRUCT_GCG_H_

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_relax.h"
#include "scip/type_sepa.h"
#include "scip/type_event.h"

#include "gcg/type_gcg.h"
#include "gcg/type_gcgpricer.h"

#ifdef __cplusplus
extern "C" {
#endif

/** GCG data structure */
struct Gcg
{
   SCIP*                origprob;            /**< SCIP data structure of origprob */
   SCIP*                masterprob;          /**< SCIP data structure of masterprob */
   SCIP*                dwmasterprob;        /**< SCIP data structure of DW masterprob */
   SCIP*                bendersmasterprob;   /**< SCIP data structure of Benders masterprob */
   SCIP_RELAX*          relax;               /**< GCG's relaxation handler */
   GCG_PRICER*          pricer;              /**< GCG pricer */
   SCIP_SEPA*           sepaorig;            /**< orig cuts separator */
   SCIP_EVENTHDLR*      mastersepacuthdlr;   /**< event handler for managing master separator cuts */
};

#ifdef __cplusplus
}
#endif

#endif /* STRUCT_GCG_H_ */
