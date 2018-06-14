/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

/**@file   pub_gcgvar.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for GCG variables
 * @author Martin Bergner
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_PUB_GCGVAR_H__
#define GCG_PUB_GCGVAR_H__

#include "type_decomp.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns TRUE or FALSE whether variable is a pricing variable or not */
extern
SCIP_Bool GCGvarIsPricing(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   );

/** returns TRUE or FALSE whether variable is a original variable or not */
extern
SCIP_Bool GCGvarIsOriginal(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   );

/** returns TRUE or FALSE whether variable is a master variable or not */
extern
SCIP_Bool GCGvarIsMaster(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   );

/** returns TRUE or FALSE whether variable is a linking variable or not */
extern
SCIP_Bool GCGoriginalVarIsLinking(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   );

/** returns TRUE or FALSE whether variable is a directly transferred variable or not */
extern
SCIP_Bool GCGoriginalVarIsTransVar(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   );

/** returns the original var of a pricing variable */
extern
SCIP_VAR* GCGpricingVarGetOriginalVar(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   );

/** adds the original var to the pricing variable */
extern
SCIP_RETCODE GCGpricingVarAddOrigVar(
   SCIP*                 scip,               /**< SCIP variable structure */
   SCIP_VAR*             pricingvar,         /**< pricing variable */
   SCIP_VAR*             origvar             /**< original pricing variable */
   );

/** returns the pricing var of an original variable */
extern
SCIP_VAR* GCGoriginalVarGetPricingVar(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns the pricing var of an original variable */
extern
void GCGoriginalVarSetPricingVar(
   SCIP_VAR*             var,                /**< SCIP variable structure */
   SCIP_VAR*             pricingvar          /**< SCIP variable structure */
   );

/** copies the pricing variable data to a master problem variable. This is used in the Benders' decomposition mode when
 * subproblems are merged into the master problem.
 */
extern
SCIP_RETCODE GCGcopyPricingvarDataToMastervar(
   SCIP*                 scip,               /**< master SCIP data structure */
   SCIP_VAR*             pricingvar,         /**< the pricing problem variable is copied from */
   SCIP_VAR*             mastervar           /**< the master variable that the vardata is copied to */
   );

/** returns the pricing variables of an linking variable */
extern
SCIP_VAR** GCGlinkingVarGetPricingVars(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   );

/** sets the pricing var of the corresponding linking variable at the specified position */
extern
void GCGlinkingVarSetPricingVar(
   SCIP_VAR*             origvar,            /**< original variable */
   int                   pricingprobnr,      /**< number of pricing problem */
   SCIP_VAR*             var                 /**< pricing variable */
   );

/** returns the number of master variables the original variable is contained in */
extern
int GCGoriginalVarGetNMastervars(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns the master variables the original variable is contained in */
extern
SCIP_VAR** GCGoriginalVarGetMastervars(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns the fraction of master variables the original variable is contained in */
extern
SCIP_Real* GCGoriginalVarGetMastervals(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns the fraction of master variables the original variable is contained in */
extern
SCIP_Real* GCGoriginalVarGetCoefs(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns the fraction of master variables the original variable is contained in */
extern
SCIP_CONS** GCGoriginalVarGetMasterconss(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** adds a coefficient of the master variable to the coefs array for the resp. constraint */
extern
SCIP_RETCODE GCGoriginalVarAddCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to add coef */
   SCIP_Real             val,                /**< coefficent to set */
   SCIP_CONS*            cons                /**< constraint the variable is in */
   );

/** adds variable to a new block, making a linkingvariable out of it, if necessary */
extern
SCIP_RETCODE GCGoriginalVarAddBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< var that is added to a block */
   int                   newblock,           /**< the new block the variable will be in */
   int                   nblocks,            /**< total number of pricing problems */
   DEC_DECMODE           mode                /**< the decomposition mode */
   );


/** returns the linking constraints */
extern
SCIP_CONS** GCGlinkingVarGetLinkingConss(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** sets the linking constraints*/
extern
void GCGlinkingVarSetLinkingCons(
   SCIP_VAR*             var,                /**< variable data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   int                   index               /**< index of pricing problem */
   );

/** returns the blocks the linking variable is in */
extern
SCIP_RETCODE GCGlinkingVarGetBlocks(
   SCIP_VAR*             var,                /**< SCIP variable structure */
   int                   nblocks,            /**< number of blocks the linking variable is in */
   int*                  blocks              /**< array to store the blocks of the linking variable */
   );

/** returns the number of blocks the linking variable is in */
extern
int GCGlinkingVarGetNBlocks(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   );

/** returns the number of coefficients of master constraints the original variable is contained in */
extern
int GCGoriginalVarGetNCoefs(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   );

/** sets the number of master variables the original variable is contained in */
extern
void GCGoriginalVarSetNCoefs(
   SCIP_VAR*             var,                /**< SCIP variable structure */
   int                   coef                /**< number of coefficient to set */
   );

/** returns TRUE or FALSE whether a master variable is a direct copy of a linking variable or not */
extern
SCIP_Bool GCGmasterVarIsLinking(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns scip instance corresponding to master variable */
extern
SCIP* GCGmasterVarGetProb(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns whether the master variable is a ray */
extern
SCIP_Bool GCGmasterVarIsRay(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns TRUE or FALSE whether a master variable is an artificial variable */
extern
SCIP_Bool GCGmasterVarIsArtificial(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns the number of original variables the master variable is contained in */
extern
int GCGmasterVarGetNOrigvars(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns the original variables the master variable is contained in */
extern
SCIP_VAR** GCGmasterVarGetOrigvars(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns the fraction of original variables the master variable is contained in */
extern
SCIP_Real* GCGmasterVarGetOrigvals(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns the number of original variables the pricing variable is contained in */
extern
int GCGpricingVarGetNOrigvars(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns the original variables the pricing variable is contained in */
extern
SCIP_VAR** GCGpricingVarGetOrigvars(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns the block of the variable */
extern
int GCGvarGetBlock(
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** sets the block of the variable */
extern
void GCGvarSetBlock(
   SCIP_VAR*             var,                /**< variable to set block for */
   int                   block               /**< block to set */
   );

/** creates the data for all variables of the original program */
extern
SCIP_RETCODE GCGcreateOrigVarsData(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** frees the data for all variables of the original program */
extern
SCIP_RETCODE GCGfreeOrigVarsData(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates the data for a variable of the original program */
extern
SCIP_RETCODE GCGorigVarCreateData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< pointer to variable object */
   );

/** returns TRUE if the linking variable is in the block, FALSE otherwise */
extern
SCIP_Bool GCGisLinkingVarInBlock(
   SCIP_VAR*             var,                /**< variable data structure */
   int                   block               /**< pricing problem number */
   );

/** determines if the master variable is in the given block */
extern
SCIP_Bool GCGisMasterVarInBlock(
   SCIP_VAR*             mastervar,          /**< master variable */
   int                   blocknr             /**< block number to check */
   );

/** informs an original variable, that a variable in the master problem was created,
 * that contains a part of the original variable.
 * Saves this information in the original variable's data
 */
extern
SCIP_RETCODE GCGoriginalVarAddMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar,            /**< original variable */
   SCIP_VAR*             var,                /**< master variable */
   SCIP_Real             val                 /**< fraction of the original variable */
   );

/* informs an original variable, that a variable in the master problem was deleted,
 * that contains a part of the original variable.
 * Update the information in the original variable's data
 */
extern
SCIP_RETCODE GCGoriginalVarRemoveMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar,            /**< original variable */
   SCIP_VAR*             var                 /**< master variable */
   );


/** creates the corresponding pricing variable for the given original variable */
extern
SCIP_RETCODE GCGoriginalVarCreatePricingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar,            /**< original variable */
   SCIP_VAR**            var                 /**< pricing variable */
   );

/** creates the corresponding pricing variable for the given original variable */
extern
SCIP_RETCODE GCGlinkingVarCreatePricingVar(
   SCIP*                 pricingscip,        /**< pricing problem SCIP data structure */
   int                   pricingprobnr,      /**< number of the pricing problem */
   SCIP_VAR*             origvar,            /**< original variable */
   SCIP_VAR**            var                 /**< pointer to store new pricing variable */
   );

/** creates the corresponding constraint in the master problem for the linking variable */
extern
SCIP_RETCODE GCGlinkingVarCreateMasterCons(
   SCIP*                 masterscip,         /**< master problem SCIP data structure */
   int                   pricingprobnr,      /**< number of the pricing problem */
   SCIP_VAR*             origvar,            /**< original variable */
   SCIP_CONS**           linkcons            /**< constraint linking pricing variables */
   );

/** creates the master var and initializes the vardata */
extern
SCIP_RETCODE GCGcreateMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origscip,           /**< original SCIP data structure */
   SCIP*                 pricingscip,        /**< pricing problem SCIP data structure */
   SCIP_VAR**            newvar,             /**< pointer to store new master variable */
   const char*           varname,            /**< new variable name */
   SCIP_Real             objcoeff,           /**< new objective coeffient */
   SCIP_VARTYPE          vartype,            /**< new variable type */
   SCIP_Bool             solisray,           /**< indicates whether new variable is a ray */
   int                   prob,               /**< number of pricing problem that created this variable */
   int                   nsolvars,           /**< number of variables in the solution */
   SCIP_Real*            solvals,            /**< values of variables in the solution */
   SCIP_VAR**            solvars,            /**< variables with non zero coefficient in the solution */
   SCIP_Bool             auxiliaryvar        /**< is new variable an Benders' auxiliary variables? */
   );

/** creates initial master variables and the vardata */
SCIP_RETCODE GCGcreateInitialMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< original variable */
   SCIP_VAR**            newvar              /**< pointer to store new variable */
   );

/** creates artificial variable and the vardata */
SCIP_RETCODE GCGcreateArtificialVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            newvar,              /**< pointer to store new variable */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             objcoef             /**< objective coefficient of artificial variable */
   );

/* adds the vardata to the auxiliary variable */
extern
SCIP_RETCODE GCGaddDataAuxiliaryVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             auxiliaryvar,       /**< the auxiliary variable */
   int                   probnumber          /**< the subproblem number */
   );

/** sets the creation node of this var */
extern
void GCGsetCreationNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< created variable */
   SCIP_Longint          creationNode        /**< node */
   );

/** returns the creation node of this var */
extern
long long int GCGgetCreationNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< created variable */
   );

/** sets the creation time of this var */
extern
void GCGsetCreationTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< created variable */
   SCIP_Real             time                /**< creation time */
   );

/** returns the creation time of this var */
extern
SCIP_Real GCGgetCreationTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< created variable */
   );

/** store pricing reduced cost call */
extern
void GCGsetRootRedcostCall(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable data structure */
   SCIP_Longint          rootredcostcall     /**< iteration at which the variable is created */
   );

/** return stored pricing reduced cost call */
extern
SCIP_Longint GCGgetRootRedcostCall(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** returns the iteration when the var was created */
extern
SCIP_Longint GCGgetIteration(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< created variable */
   );

/** sets the iteration when the var was created */
extern
void GCGsetIteration(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< created variable */
   SCIP_Longint          iteration           /**< iteration that this var was created */
   );

/** store gap */
void GCGsetGap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable data structure */
   SCIP_Real             gap                 /**< present gap when variable is created */
   );

/** return stored gap */
SCIP_Real GCGgetGap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** store reduced cost */
void GCGsetRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable data structure */
   SCIP_Real             redcost             /**< reduced cost of the variable at creation */
   );

/** return stored reduced cost */
SCIP_Real GCGgetRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable data structure */
   );

/** updates the statistics part of the variable */
void GCGupdateVarStatistics(
    SCIP*                scip,               /**< master SCIP data structure */
    SCIP*                origprob,           /**< original SCIP data structure */
    SCIP_VAR*            newvar,             /**< new variable for statistic update */
    SCIP_Real            redcost             /**< reduced cost of the variable */
   );

/** prints the given variable: name, type (original, master or pricing) block number,
 * and the list of all variables related to the given variable */
extern
void GCGprintVar(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File to write information to, or NULL for stdout */
   SCIP_VAR*             var                 /**< variable that should be printed */
   );

#ifdef __cplusplus
}
#endif

#endif
