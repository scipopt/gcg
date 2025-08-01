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

/**@file    struct_extendedmasterconsdata.h
 * @ingroup DATASTRUCTURES
 * @brief   data structures for GCG extended master cons data
 * @author  Til Mohr
 * @author  Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_EXTENDEDMASTERCONSDATA_H_
#define GCG_STRUCT_EXTENDEDMASTERCONSDATA_H_

#include <scip/def.h>
#include <scip/type_cons.h>
#include <scip/type_lp.h>
#include <scip/type_var.h>
#include "gcg/type_extendedmasterconsdata.h"
#include "gcg/type_branchgcg.h"
#include "gcg/type_mastersepacut.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data for a pricing problem modification */
struct GCG_PricingModification
{
   int                        blocknr;            /**< block number of the extended master cons */
   SCIP_VAR*                  coefvar;            /**< variable in the pricing problem inferred from the extended master cons
                                                    * always has the objective coefficient of the negated dual value of the extended master cons
                                                    * its solution value corresponds to the coefficient of the new mastervariable in the extended master cons */
   SCIP_VAR**                 additionalvars;     /**< array of additional variables with no objective coefficient in the pricing programs inferred from the extended master cons */
   int                        nadditionalvars;    /**< number of additional variables in the pricing programs */
   SCIP_CONS**                additionalconss;    /**< array of additional constraints in the pricing programs inferred from the extended master cons */
   int                        nadditionalconss;   /**< number of additional constraints in the pricing programs */
};

/** cut of the extended master cons */
union GCG_ExtendedMasterCons
{
   SCIP_CONS*                 cons;                /**< constraint in the master problem that represents the extended master cons, iff type == Cons */
   SCIP_ROW*                  row;                 /**< row in the master problem that represents the extended master cons, iff type == Row */
};

/** data for extended master conss */
struct GCG_ExtendedMasterConsData
{
   GCG_EXTENDEDMASTERCONSTYPE       type;                      /**< type of the extended master cons */
   GCG_EXTENDEDMASTERCONS           cons;                      /**< constraint or row in the master problem that represents the extended master cons */
   GCG_PRICINGMODIFICATION**        pricingmodifications;      /**< array of pricing modifications for the extended master cons */
   int                              npricingmodifications;     /**< number of pricing modifications for the extended master cons */
   union {
      GCG_BRANCHCONSDATA*           branchconsdata;            /**< branchconsdata in case the extended master cons is of type GCG_EXTENDEDMASTERCONSTYPE_BRANCH_CONS */
      GCG_SEPARATORMASTERCUT*       sepamastercut;             /**< sepamastercut in case the extended master cons is of type GCG_EXTENDEDMASTERCONSTYPE_SEPA_ROW */
   } data;                                                     /**< data required to calculate the coefficient of a column solution */
};

#ifdef __cplusplus
}
#endif

#endif
