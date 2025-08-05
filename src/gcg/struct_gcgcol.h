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

/**@file   struct_gcgcol.h
 * @ingroup DATASTRUCTURES
 * @brief  data structure to store columns (solutions from a pricing problem)
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_GCGCOL_H_
#define GCG_STRUCT_GCGCOL_H_

#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/scip.h"

#include "gcg/type_gcgcol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Column data structure */
struct GCG_Col
{
   SCIP*                pricingprob;        /**< SCIP data structure (pricing problem)*/
   int                  probnr;             /**< number of corresponding pricing problem */
   SCIP_VAR**           vars;               /**< (sorted) array of variables of corresponding pricing problem */
   SCIP_Real*           vals;               /**< array of solution values (belonging to vars) */
   int                  nvars;              /**< number of variables */
   int                  maxvars;            /**< capacity of vars */
   SCIP_VAR**           inferredpricingvars;/**< (sorted) array of inferrred (coefficient) pricing variables */
   SCIP_Real*           inferredpricingvals;/**< inferred (coefficient) pricing variables solution values */
   int                  ninferredpricingvars;/**< number of inferred pricing variables */
   int                  maxinferredpricingvars;/**< capacity of inferred pricing variables array */
   SCIP_Bool            isray;              /**< is the column a ray? */
   SCIP_Real            redcost;            /**< last known reduced cost */
   int                  age;                /**< age of column (number of iterations since it was created;
                                                 each time reduced cost are calculated counts as an interation) */
   int                  pos;                /**< position in pricestore or column pool (or -1) */
   SCIP_Real*           mastercoefs;        /**< array of master coefficients */
   int                  nmastercoefs;       /**< number of master coefficients */
   int                  maxmastercoefs;     /**< capacity of mastercoefs */
   SCIP_Real*           originalsepamastercuts;/**< array of original seperator cut coefficients in the master problem */
   int                  noriginalsepamastercuts;/**< number of original seperator cut coefficients in the master problem */
   int                  maxoriginalsepamastercuts;/**< capacity of originalsepamastercuts */
   SCIP_Real*           mastersepacutcoeffs;        /**< arrays of master separator cut coefficients */
   int                  nmastersepacutcoeffs;       /**< number of master separator cut coefficients */
   int                  mastersepacutscoeffssize;   /**< size of array of master separator cut coefficients */
   SCIP_Real            norm;               /**< norm of the coefficients in the master */
   int*                 linkvars;           /**< array of indices of variables in var-array which are linking variables */
   int                  nlinkvars;          /**< number of variables in var-array which are linking variables */
   int                  maxlinkvars;        /**< capacity of linkvars */
   SCIP_Bool            initcoefs;          /**< returns if mastercoefs and linkvars have been computed */
};

#ifdef __cplusplus
}
#endif

#endif /* STRUCT_GCGCOL_H_ */
