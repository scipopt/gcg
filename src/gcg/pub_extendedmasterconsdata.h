/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
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
#include "gcg/def.h"
#include "scip/scip.h"

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
   SCIP*                         scip,                         /**< SCIP data structure */
   GCG_PRICINGMODIFICATION*      pricingmodification,          /**< pointer to store the pricing modification */
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
   SCIP*                         scip,                         /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA**  extendedmasterconsdata,       /**< pointer to store the extended master cons data */
   SCIP_CONS*                    cons,                         /**< constraint in the master problem that represents the extended master cons */
   GCG_PRICINGMODIFICATION*      pricingmodifications,         /**< pricing modifications for the extended master cons */
   int                           npricingmodifications,        /**< number of pricing modifications for the extended master cons */
   void*                         data,                         /**< any data that might be required to calculate the coefficient of a column solution */
   GCG_DECL_EXTENDEDMASTERCONSGETCOEFF((*extendedmasterconsGetCoeff))   /**< callback to calculate the coefficient of a column solution */
   );

/** create an extended master cons, taking ownership over pricingmodifications */
GCG_EXPORT
SCIP_RETCODE GCGextendedmasterconsCreateFromRow(
   SCIP*                         scip,                         /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA**  extendedmasterconsdata,       /**< pointer to store the extended master cons data */
   SCIP_ROW*                     row,                          /**< row in the master problem that represents the extended master cons */
   GCG_PRICINGMODIFICATION*      pricingmodifications,         /**< pricing modifications for the extended master cons */
   int                           npricingmodifications,        /**< number of pricing modifications for the extended master cons */
   void*                         data,                         /**< any data that might be required to calculate the coefficient of a column solution */
   GCG_DECL_EXTENDEDMASTERCONSGETCOEFF((*extendedmasterconsGetCoeff))   /**< callback to calculate the coefficient of a column solution */
   );

/** free an extended master cons */
GCG_EXPORT
SCIP_RETCODE GCGextendedmasterconsFree(
   SCIP*                         scip,                         /**< SCIP data structure */
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
   SCIP*                         masterscip,                   /**< master scip */
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
   SCIP*                         masterscip,                   /**< master scip */
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata,       /**< extended master cons data */
   int                           blocknr                       /**< block number */
   );

#ifdef NDEBUG
#define GCGextendedmasterconsGetPricingModifications(extendedmasterconsdata) extendedmasterconsdata->pricingmodifications
#else
/** get the pricing modifications for the extended master cons */
GCG_EXPORT
GCG_PRICINGMODIFICATION* GCGextendedmasterconsGetPricingModifications(
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
   GCG_PRICINGMODIFICATION       pricingmodification,          /**< pricing modification */
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
   SCIP*                         scip,                         /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*   extendedmasterconsdata        /**< extended master cons data */
   );

/** get the rhs of the extended master cons */
GCG_EXPORT
SCIP_Real GCGextendedmasterconsGetRhs(
   SCIP*                  scip,                                /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata      /**< extended master cons data */
   );

/** get the constant of the extended master cons (always returns 0 if extended master cons is a constraint, returns constant of row otherwise) */
GCG_EXPORT
SCIP_Real GCGextendedmasterconsGetConstant(
   SCIP*                            scip,                      /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata     /**< extended master cons data */
   );

/** get number of nonzero entries in the extended master cons */
GCG_EXPORT
int GCGextendedmasterconsGetNNonz(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*     extendedmasterconsdata       /**< extended master cons data */
   );

/** get array of columns with nonzero entries */
GCG_EXPORT
SCIP_COL** GCGextendedmasterconsGetCols(
   SCIP*                            scip,                         /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata        /**< extended master cons data */
   );

/** get array of coefficients with nonzero entries */
GCG_EXPORT
SCIP_Real* GCGextendedmasterconsGetVals(
   SCIP*                            scip,                         /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata        /**< extended master cons data */
   );

#ifdef NDEBUG
#define GCGextendedmasterconsGetData(extendedmasterconsdata) extendedmasterconsdata->data
#else
/** get the additional data */
GCG_EXPORT
void* GCGextendedmasterconsGetData(
   GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata        /**< extended master cons data */
   );
#endif

/** calculate the coefficient of a column solution in the extended master cons */
GCG_EXPORT
SCIP_Real GCGextendedmasterconsGetCoeff(
   SCIP*                            scip,                         /**< SCIP data structure */
   GCG_EXTENDEDMASTERCONSDATA*      extendedmasterconsdata,       /**< extended master cons data */
   SCIP_VAR**                       solvars,                      /**< array of column solution variables */
   SCIP_Real*                       solvals,                      /**< array of column solution values */
   int                              nsolvars,                     /**< number of column solution variables and values */
   int                              probnr                        /**< the pricing problem that the column belongs to */
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
