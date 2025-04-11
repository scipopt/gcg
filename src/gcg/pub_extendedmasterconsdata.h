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

/**@file    pub_extendedmasterconsdata.h
 * @ingroup PUBLICCOREAPI
 * @brief   public methods for interacting with GCG_EXTENDEDMASTERCONSDATA
 * @author  Til Mohr
 * @author  Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_PUB_EXTENDEDMASTERCONSDATA_H_
#define GCG_PUB_EXTENDEDMASTERCONSDATA_H_

#include "gcg/type_extendedmasterconsdata.h"

#include "scip/scip.h"
#include "gcg/def.h"
#include "gcg/type_gcg.h"

#ifdef NDEBUG
#include "gcg/struct_extendedmasterconsdata.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @ingroup GCG_EXTENDEDMASTERCONSDATA
 *
 * @{
 */

/** create a pricing modification, taking ownership over additionalvars and additionalcons */
GCG_EXPORT
SCIP_RETCODE GCGpricingmodificationCreate(
   GCG*                          gcg,                          /**< GCG data structure */
   GCG_PRICINGMODIFICATION**     pricingmodification,          /**< pointer to store the pricing modification */
   int                           blocknr,                      /**< block number of the extended master cons */
   SCIP_VAR*                     coefvar,                      /**< variable in the pricing problem inferred from the extended master cons
                                                                * always has the objective coefficient of the negated dual value of the extended master cons
                                                                * its solution value corresponds to the coefficient of the new mastervariable in the extended master cons */
   SCIP_VAR**                    additionalvars,               /**< array of additional variables with no objective coefficient in the pricing programs inferred from the extended master cons */
   int                           nadditionalvars,              /**< number of additional variables in the pricing programs */
   SCIP_CONS**                   additionalconss,              /**< array of additional constraints in the pricing programs inferred from the extended master cons */
   int                           nadditionalconss              /**< number of additional constraints in the pricing programs */
   );

/** create an extended master cons, taking ownership over pricingmodifications */
GCG_EXPORT
SCIP_RETCODE GCGextendedmasterconsCreateFromCons(
   GCG*                             gcg,                          /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA**     extendedmasterconsdata,       /**< pointer to store the extended master cons data */
   GCG_EXTENDEDMASTERCONSTYPE       type,                         /**< type of the extended master cons */
   SCIP_CONS*                       cons,                         /**< constraint in the master problem that represents the extended master cons */
   GCG_PRICINGMODIFICATION**        pricingmodifications,         /**< pricing modifications for the extended master cons */
   int                              npricingmodifications,        /**< number of pricing modifications for the extended master cons */
   GCG_EXTENDEDMASTERCONSDATADATA*  data                          /**< any data that might be required to calculate the coefficient of a column solution */
   );

/** create an extended master cons, taking ownership over pricingmodifications */
GCG_EXPORT
SCIP_RETCODE GCGextendedmasterconsCreateFromRow(
   GCG*                             gcg,                          /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA**     extendedmasterconsdata,       /**< pointer to store the extended master cons data */
   GCG_EXTENDEDMASTERCONSTYPE       type,                         /**< type of the extended master cons */
   SCIP_ROW*                        row,                          /**< row in the master problem that represents the extended master cons */
   GCG_PRICINGMODIFICATION**        pricingmodifications,         /**< pricing modifications for the extended master cons */
   int                              npricingmodifications,        /**< number of pricing modifications for the extended master cons */
   GCG_EXTENDEDMASTERCONSDATADATA*  data                          /**< any data that might be required to calculate the coefficient of a column solution */
   );

/** free an extended master cons */
GCG_EXPORT
SCIP_RETCODE GCGextendedmasterconsFree(
   GCG*                          gcg,                          /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA**  extendedmasterconsdata        /**< pointer to the extended master cons data */
   );

/** determine whether the extendedmasterconsdata is active in the masterscip */
GCG_EXPORT
SCIP_Bool GCGextendedmasterconsIsActive(
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata        /**< extended master cons data */
   );

/** add a new variable along with its coefficient to the extended master cons */
GCG_EXPORT
SCIP_RETCODE GCGextendedmasterconsAddMasterVar(
   GCG*                          gcg,                          /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata,       /**< extended master cons data */
   SCIP_VAR*                     var,                          /**< variable to add */
   SCIP_Real                     coef                          /**< coefficient of the variable */
   );

#ifdef NDEBUG
#define GCGextendedmasterconsGetCons(extendedmasterconsdata) extendedmasterconsdata->cons.cons
#else
/** get the constraint that is the extended master cons */
GCG_EXPORT
SCIP_CONS* GCGextendedmasterconsGetCons(
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata        /**< extended master cons data */
   );
#endif

#ifdef NDEBUG
#define GCGextendedmasterconsGetRow(extendedmasterconsdata) extendedmasterconsdata->cons.row
#else
/** get the row that is the extended master cons */
GCG_EXPORT
SCIP_ROW* GCGextendedmasterconsGetRow(
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata        /**< extended master cons data */
   );
#endif

#ifdef NDEBUG
#define GCGpricingmodificationGetCoefVar(pricingmodification) pricingmodification->coefvar
#else
/** get the variable that determines the coefficient of a column in the extended master cons */
GCG_EXPORT
SCIP_VAR* GCGpricingmodificationGetCoefVar(
   GCG_PRICINGMODIFICATION*      pricingmodification           /**< pricing modification */
   );
#endif

#ifdef NDEBUG
#define GCGpricingmodificationGetAdditionalVars(pricingmodification) pricingmodification->additionalvars
#else
/** get the additional variables that are inferred by the extended master cons */
GCG_EXPORT
SCIP_VAR** GCGpricingmodificationGetAdditionalVars(
   GCG_PRICINGMODIFICATION*      pricingmodification           /**< pricing modification */
   );
#endif

#ifdef NDEBUG
#define GCGpricingmodificationGetNAdditionalVars(pricingmodification) pricingmodification->nadditionalvars
#else
/** get the number of additional variables that are inferred by the extended master cons */
GCG_EXPORT
int GCGpricingmodificationGetNAdditionalVars(
   GCG_PRICINGMODIFICATION*      pricingmodification           /**< pricing modification */
   );
#endif

#ifdef NDEBUG
#define GCGpricingmodificationGetAdditionalConss(pricingmodification) pricingmodification->additionalconss
#else
/** get the additional constraints that are inferred by the extended master cons */
GCG_EXPORT
SCIP_CONS** GCGpricingmodificationGetAdditionalConss(
   GCG_PRICINGMODIFICATION*      pricingmodification           /**< pricing modification */
   );
#endif

#ifdef NDEBUG
#define GCGpricingmodificationGetNAdditionalConss(pricingmodification) pricingmodification->nadditionalconss
#else
/** get the number of additional constraints that are inferred by the extended master cons */
GCG_EXPORT
int GCGpricingmodificationGetNAdditionalConss(
   GCG_PRICINGMODIFICATION*      pricingmodification           /**< pricing modification */
   );
#endif

/** get the pricing modification for a block, if exists, else NULL */
GCG_EXPORT
GCG_PRICINGMODIFICATION* GCGextendedmasterconsGetPricingModification(
   GCG*                          gcg,                          /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata,       /**< extended master cons data */
   int                           blocknr                       /**< block number */
   );

#ifdef NDEBUG
#define GCGextendedmasterconsGetPricingModifications(extendedmasterconsdata) extendedmasterconsdata->pricingmodifications
#else
/** get the pricing modifications for the extended master cons */
GCG_EXPORT
GCG_PRICINGMODIFICATION** GCGextendedmasterconsGetPricingModifications(
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata        /**< extended master cons data */
   );
#endif

#ifdef NDEBUG
#define GCGextendedmasterconsGetNPricingModifications(extendedmasterconsdata) extendedmasterconsdata->npricingmodifications
#else
/** get the number of pricing modifications for the extended master cons */
GCG_EXPORT
int GCGextendedmasterconsGetNPricingModifications(
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata        /**< extended master cons data */
   );
#endif

/** check whether a given variable is a coefficient variable of a given pricing modification */
GCG_EXPORT
SCIP_Bool GCGpricingmodificationIsCoefVar(
   GCG_PRICINGMODIFICATION*      pricingmodification,          /**< pricing modification */
   SCIP_VAR*                     var                           /**< variable to check */
   );

/** check whether a given variable is a coefficient variable of a given extended master cons */
GCG_EXPORT
SCIP_Bool GCGextendedmasterconsIsCoefVar(
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata,       /**< extended master cons data */
   SCIP_VAR*                     var                           /**< variable to check */
   );

/** get name of the extended master cons */
GCG_EXPORT
const char* GCGextendedmasterconsGetName(
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata        /**< extended master cons data */
   );

/** get the lhs of the extended master cons */
GCG_EXPORT
SCIP_Real GCGextendedmasterconsGetLhs(
   GCG*                          gcg,                          /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata        /**< extended master cons data */
   );

/** get the rhs of the extended master cons */
GCG_EXPORT
SCIP_Real GCGextendedmasterconsGetRhs(
   GCG*                            gcg,                        /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata      /**< extended master cons data */
   );

/** get the constant of the extended master cons (always returns 0 if extended master cons is a constraint, returns constant of row otherwise) */
GCG_EXPORT
SCIP_Real GCGextendedmasterconsGetConstant(
   GCG*                             gcg,                       /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata     /**< extended master cons data */
   );

/** get number of nonzero entries in the extended master cons */
GCG_EXPORT
int GCGextendedmasterconsGetNNonz(
   GCG*                            gcg,                         /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   );

/** get array of columns with nonzero entries */
GCG_EXPORT
SCIP_COL** GCGextendedmasterconsGetCols(
   GCG*                             gcg,                          /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata        /**< extended master cons data */
   );

/** get array of coefficients with nonzero entries */
GCG_EXPORT
SCIP_Real* GCGextendedmasterconsGetVals(
   GCG*                             gcg,                          /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata        /**< extended master cons data */
   );

#ifdef NDEBUG
#define GCGextendedmasterconsGetData(extendedmasterconsdata) extendedmasterconsdata->data
#else
/** get the additional data */
GCG_EXPORT
GCG_EXTENDEDMASTERCONSDATADATA* GCGextendedmasterconsGetData(
   GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata        /**< extended master cons data */
   );
#endif

/** calculate the coefficient of a column solution in the extended master cons */
GCG_EXPORT
SCIP_RETCODE GCGextendedmasterconsGetCoeff(
   GCG*                             gcg,                          /**< GCG data structure */
   GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata,       /**< extended master cons data */
   SCIP_VAR**                       solvars,                      /**< array of column solution variables */
   SCIP_Real*                       solvals,                      /**< array of column solution values */
   int                              nsolvars,                     /**< number of column solution variables and values */
   int                              probnr,                       /**< the pricing problem that the column belongs to */
   SCIP_Real*                       coeff                         /**< pointer to store the coefficient */
   );

#ifdef NDEBUG
#define GCGextendedmasterconsGetType(extendedmasterconsdata) extendedmasterconsdata->type
#else
/** gets the type of the extended master cons */
GCG_EXPORT
GCG_EXTENDEDMASTERCONSTYPE GCGextendedmasterconsGetType(
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata           /**< extended master cons data */
   );
#endif

/**@} */
#ifdef __cplusplus
}
#endif

#endif
