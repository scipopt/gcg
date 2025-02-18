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

/**@file    pub_mastercutdata.h
 * @ingroup PUBLICCOREAPI
 * @brief   public methods for interacting with GCG_MASTERCUTDATA
 * @author  Til Mohr
 * @author  Erik Muehmer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_PUB_MASTERCUTDATA_H_
#define GCG_PUB_MASTERCUTDATA_H_

#include "struct_mastercutdata.h"
#include "type_mastercutdata.h"
#include "def.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @ingroup GCG_MASTERCUTDATA
 *
 * @{
 */

/** create a pricing modification, taking ownership over additionalvars and additionalcons */
GCG_EXPORT
SCIP_RETCODE GCGpricingmodificationCreate(
   SCIP*                  scip,                 /**< SCIP data structure */
   GCG_PRICINGMODIFICATION* pricingmodification,/**< pointer to store the pricing modification */
   int                    blocknr,              /**< block number of the master cut */
   SCIP_VAR*              coefvar,              /**< variable in the pricing problem inferred from the master cut
                                                  * always has the objective coefficient of the negated dual value of the master cut
                                                  * its solution value corresponds to the coefficient of the new mastervariable in the master cut */
   SCIP_VAR**             additionalvars,       /**< array of additional variables with no objective coefficient in the pricing programs inferred from the master cut */
   int                    nadditionalvars,      /**< number of additional variables in the pricing programs */
   SCIP_CONS**            additionalconss,      /**< array of additional constraints in the pricing programs inferred from the master cut */
   int                    nadditionalconss      /**< number of additional constraints in the pricing programs */
   );

/** create a master cut, taking ownership over pricingmodifications */
GCG_EXPORT
SCIP_RETCODE GCGmastercutCreateFromCons(
   SCIP*                  scip,                 /**< SCIP data structure */
   GCG_MASTERCUTDATA**    mastercutdata,        /**< pointer to store the mastercut data */
   SCIP_CONS*             cons,                 /**< constraint in the master problem that represents the master cut */
   GCG_PRICINGMODIFICATION* pricingmodifications,/**< pricing modifications for the master cut */
   int                    npricingmodifications,/**< number of pricing modifications for the master cut */
   void*                  data,                 /**< any data that might be required to calculate the coefficient of a column solution */
   GCG_DECL_MASTERCUTGETCOEFF((*mastercutGetCoeff))/**< callback to calculate the coefficient of a column solution */
   );

/** create a master cut, taking ownership over pricingmodifications */
GCG_EXPORT
SCIP_RETCODE GCGmastercutCreateFromRow(
   SCIP*                  scip,                 /**< SCIP data structure */
   GCG_MASTERCUTDATA**    mastercutdata,        /**< pointer to store the mastercut data */
   SCIP_ROW*              row,                  /**< row in the master problem that represents the master cut */
   GCG_PRICINGMODIFICATION* pricingmodifications,/**< pricing modifications for the master cut */
   int                    npricingmodifications,/**< number of pricing modifications for the master cut */
   void*                  data,                 /**< any data that might be required to calculate the coefficient of a column solution */
   GCG_DECL_MASTERCUTGETCOEFF((*mastercutGetCoeff))/**< callback to calculate the coefficient of a column solution */
   );

/** free a master cut */
GCG_EXPORT
SCIP_RETCODE GCGmastercutFree(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA**    mastercutdata        /**< pointer to the mastercut data */
   );

/** determine whether the mastercutdata is active in the masterscip */
GCG_EXPORT
SCIP_Bool GCGmastercutIsActive(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );

/** add a new variable along with its coefficient to the mastercut */
GCG_EXPORT
SCIP_RETCODE GCGmastercutAddMasterVar(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   SCIP_VAR*              var,                /**< variable to add */
   SCIP_Real              coef                /**< coefficient of the variable */
   );

#ifdef NDEBUG
#define GCGmastercutGetCons(mastercutdata) mastercutdata->cut.cons
#else
/** get the constraint that is the master cut */
GCG_EXPORT
SCIP_CONS* GCGmastercutGetCons(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );
#endif

#ifdef NDEBUG
#define GCGmastercutGetRow(mastercutdata) mastercutdata->cut.row
#else
/** get the row that is the master cut */
GCG_EXPORT
SCIP_ROW* GCGmastercutGetRow(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );
#endif

/** get the variable that determines the coefficient of a column in the master cut */
GCG_EXPORT
SCIP_VAR* GCGpricingmodificationGetCoefVar(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   );

/** get the additional variables that are inferred by the master cut */
GCG_EXPORT
SCIP_VAR** GCGpricingmodificationGetAdditionalVars(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   );

/** get the number of additional variables that are inferred by the master cut */
GCG_EXPORT
int GCGpricingmodificationGetNAdditionalVars(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   );

/** get the additional constraints that are inferred by the master cut */
GCG_EXPORT
SCIP_CONS** GCGpricingmodificationGetAdditionalConss(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   );

/** get the number of additional constraints that are inferred by the master cut */
GCG_EXPORT
int GCGpricingmodificationGetNAdditionalConss(
   GCG_PRICINGMODIFICATION* pricingmodification /**< pricing modification */
   );

/** get the pricing modification for a block, if exists, else NULL */
GCG_EXPORT
GCG_PRICINGMODIFICATION* GCGmastercutGetPricingModification(
   SCIP*                  masterscip,         /**< master scip */
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   int                    blocknr             /**< block number */
   );

#ifdef NDEBUG
#define GCGmastercutGetPricingModifications(mastercutdata) mastercutdata->pricingmodifications
#else
/** get the pricing modifications for the master cut */
GCG_EXPORT
GCG_PRICINGMODIFICATION* GCGmastercutGetPricingModifications(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );
#endif

#ifdef NDEBUG
#define GCGmastercutGetNPricingModifications(mastercutdata) mastercutdata->npricingmodifications
#else
/** get the number of pricing modifications for the master cut */
GCG_EXPORT
int GCGmastercutGetNPricingModifications(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );
#endif

/** check whether a given variable is a coefficient variable of a given pricing modification */
GCG_EXPORT
SCIP_Bool GCGpricingmodificationIsCoefVar(
   GCG_PRICINGMODIFICATION pricingmodification, /**< pricing modification */
   SCIP_VAR*              var                 /**< variable to check */
   );

/** check whether a given variable is a coefficient variable of a given mastercut */
GCG_EXPORT
SCIP_Bool GCGmastercutIsCoefVar(
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   SCIP_VAR*              var                 /**< variable to check */
   );

/** get name of the mastercut */
GCG_EXPORT
const char* GCGmastercutGetName(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );

/** get the lhs of the mastercut */
GCG_EXPORT
SCIP_Real GCGmastercutGetLhs(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );

/** get the rhs of the mastercut */
GCG_EXPORT
SCIP_Real GCGmastercutGetRhs(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );

/** get the constant of the mastercut (always returns 0 if mastercut is a constraint, returns constant of row otherwise) */
GCG_EXPORT
SCIP_Real GCGmastercutGetConstant(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );

/** get number of nonzero entries in the mastercut */
GCG_EXPORT
int GCGmastercutGetNNonz(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );

/** get array of columns with nonzero entries */
GCG_EXPORT
SCIP_COL** GCGmastercutGetCols(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );

/** get array of coefficients with nonzero entries */
GCG_EXPORT
SCIP_Real* GCGmastercutGetVals(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );

/** get the additional data */
GCG_EXPORT
void* GCGmastercutGetData(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );

/** calculate the coefficient of a column solution in the master cut */
GCG_EXPORT
SCIP_Real GCGmastercutGetCoeff(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_MASTERCUTDATA*     mastercutdata,      /**< mastercut data */
   SCIP_VAR**             solvars,            /**< array of column solution variables */
   SCIP_Real*             solvals,            /**< array of column solution values */
   int                    nsolvars,           /**< number of column solution variables and values */
   int                    probnr              /**< the pricing problem that the column belongs to */
   );

#ifdef NDEBUG
#define GCGmastercutGetType(mastercutdata) mastercutdata->type
#else
/** gets the type of the master cut */
GCG_EXPORT
GCG_MASTERCUTTYPE GCGmastercutGetType(
   GCG_MASTERCUTDATA*     mastercutdata       /**< mastercut data */
   );
#endif

/**@} */
#ifdef __cplusplus
}
#endif

#endif
