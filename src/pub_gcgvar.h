/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dec_borderheur.h
 * @brief  borderheur presolver
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GCG_PUB_GCGVAR_H__
#define __GCG_PUB_GCGVAR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Returns TRUE or FALSE whether variable is a pricing variable or not */
extern
SCIP_Bool GCGvarIsPricing(
   SCIP_VAR* var          /**< SCIP variable */
   );

/** Returns TRUE or FALSE whether variable is a original variable or not */
extern
SCIP_Bool GCGvarIsOriginal(
   SCIP_VAR* var          /**< SCIP variable structure */
   );

/** Returns TRUE or FALSE whether variable is a master variable or not */
extern
SCIP_Bool GCGvarIsMaster(
   SCIP_VAR* var          /**< SCIP variable structure */
   );

/** Returns TRUE or FALSE whether variable is a linking variable or not */
extern
SCIP_Bool GCGvarIsLinking(
   SCIP_VAR* var /**< SCIP variable structure */
   );

/** Returns the original var of a pricing variable */
extern
SCIP_VAR* GCGpricingVarGetOriginalVar(
   SCIP_VAR* var /**< SCIP variable structure */
   );

/** Adds the original var to the pricing variable */
extern
SCIP_RETCODE GCGpricingVarAddOrigVar(
   SCIP* scip, /**< SCIP variable structure */
   SCIP_VAR* pricingvar,
   SCIP_VAR* origvar
   );

/** Returns the pricing var of an original variable */
extern
SCIP_VAR* GCGoriginalVarGetPricingVar(
   SCIP_VAR* var /**< SCIP variable structure */
   );

/** Returns the pricing var of an original variable */
extern
void GCGoriginalVarSetPricingVar(
   SCIP_VAR* var, /**< SCIP variable structure */
   SCIP_VAR* pricingvar /**< SCIP variable structure */
   );

/** Returns the pricing variables of an linking variable */
extern
SCIP_VAR** GCGlinkingVarGetPricingVars(
   SCIP_VAR* var /**< SCIP variable structure */
   );

/** sets the pricing var of the corresponding linking variable at the specified position */
extern
void GCGlinkingVarSetPricingVar(
   SCIP_VAR* origvar,
   int pricingprobnr,
   SCIP_VAR* var
   );

/** Returns the number of master variables the original variable is contained in */
extern
int GCGoriginalVarGetNMastervars(
   SCIP_VAR* var
   );

/** Returns the master variables the original variable is contained in */
extern
SCIP_VAR** GCGoriginalVarGetMastervars(
   SCIP_VAR* var
   );

/** Returns the fraction of master variables the original variable is contained in */
extern
SCIP_Real* GCGoriginalVarGetMastervals(
   SCIP_VAR* var
   );

/** Returns the fraction of master variables the original variable is contained in */
extern
SCIP_Real* GCGoriginalVarGetCoefs(
   SCIP_VAR* var
   );

/** Returns the fraction of master variables the original variable is contained in */
extern
SCIP_CONS** GCGoriginalVarGetMasterconss(
   SCIP_VAR* var
   );

/** adds a coefficient of the master variable to the coefs array for the resp. constraint */
extern
SCIP_RETCODE GCGoriginalVarAddCoef(
   SCIP* scip,         /**< SCIP data structure */
   SCIP_VAR* var,      /**< variable to add coef */
   SCIP_Real val,      /**< coefficent to set */
   SCIP_CONS* cons     /**< constraint the variable is in */
   );

/** adds variable to a new block, making a linkingvariable out of it, if necessary */
extern
SCIP_RETCODE GCGoriginalVarAddBlock(
   SCIP* scip,    /**< SCIP data structure */
   SCIP_VAR* var, /**< var that is added to a block */
   int newblock   /**< The new block the variable will be in */
   );


/** Returns the fraction of master variables the original variable is contained in */
extern
SCIP_CONS** GCGlinkingVarGetLinkingConss(
   SCIP_VAR* var
   );

/** Returns the fraction of master variables the original variable is contained in */
extern
void GCGlinkingVarSetLinkingCons(
   SCIP_VAR* var,
   SCIP_CONS* cons,
   int i
   );

/** returns the number of blocks the linking variable is in */
extern
int GCGlinkingVarGetNBlocks(
   SCIP_VAR* var /**< SCIP variable structure */
   );

/** Returns the number of master variables the original variable is contained in */
extern
int GCGoriginalVarGetNCoefs(
   SCIP_VAR* var /**< SCIP variable structure */
   );

/** Sets the number of master variables the original variable is contained in */
extern
void GCGoriginalVarSetNCoefs(
   SCIP_VAR* var, /**< SCIP variable structure      */
   int coef       /**< number of coefficient to set */
   );

/** Returns whether the  master variable is a ray */
extern
SCIP_Bool GCGmasterVarIsRay(
   SCIP_VAR* var
   );

/** Returns the number of original variables the master variable is contained in */
extern
int GCGmasterVarGetNOrigvars(
   SCIP_VAR* var
   );

/** Returns the original variables the master variable is contained in */
extern
SCIP_VAR** GCGmasterVarGetOrigvars(
   SCIP_VAR* var
   );

/** Returns the fraction of original variables the master variable is contained in */
extern
SCIP_Real* GCGmasterVarGetOrigvals(
   SCIP_VAR* var
   );

/** Returns the number of original variables the pricing variable is contained in */
extern
int GCGpricingVarGetNOrigvars(
   SCIP_VAR* var
   );

/** Returns the original variables the pricing variable is contained in */
extern
SCIP_VAR** GCGpricingVarGetOrigvars(
   SCIP_VAR* var
   );

/** Returns the block of the variable */
extern
int GCGvarGetBlock(
   SCIP_VAR* var
   );

/** Sets the block of the variable */
extern
void GCGvarSetBlock(
   SCIP_VAR* var, /**< Variable to set block for */
   int block      /**< block to set */
   );

/** creates the data for all variables of the original program */
extern
SCIP_RETCODE GCGcreateOrigVarsData(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates the data for a variable of the original program */
extern
SCIP_RETCODE GCGorigVarCreateData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< pointer to variable object */
   );

/** Returns TRUE if the linking variable is in the block, FALSE otherwise */
extern
SCIP_Bool GCGisLinkingVarInBlock(
   SCIP_VAR* var,
   int block
   );

/* informs an original variable, that a variable in the master problem was created,
 * that contains a part of the original variable.
 * Saves this information in the original variable's data
 * @todo this method needs a little love
 */
extern
SCIP_RETCODE GCGoriginalVarAddMasterVar(
   SCIP*                 scip,                  /**< SCIP data structure                */
   SCIP_VAR*             origvar,               /**< Original variable                  */
   SCIP_VAR*             var,                   /**< Master variable                    */
   SCIP_Real             val                    /**< Fraction of the original variable  */
   );

/* informs an original variable, that a variable in the master problem was deleted,
 * that contains a part of the original variable.
 * Update the information in the original variable's data
 * @todo this method needs a little love
 */
extern
SCIP_RETCODE GCGoriginalVarRemoveMasterVar(
   SCIP*                 scip,                  /**< SCIP data structure                */
   SCIP_VAR*             origvar,               /**< Original variable                  */
   SCIP_VAR*             var                    /**< Master variable                    */
   );



/** creates the corresponding pricing variable for the given original variable */
extern
SCIP_RETCODE GCGoriginalVarCreatePricingVar(
   SCIP* scip,
   SCIP_VAR* origvar,
   SCIP_VAR** var
   );

/** creates the corresponding pricing variable for the given original variable */
extern
SCIP_RETCODE GCGlinkingVarCreatePricingVar(
   SCIP* masterscip,
   SCIP* pricingscip,
   int pricingprobnr,
   SCIP_VAR* origvar,
   SCIP_VAR** var,
   SCIP_CONS** linkcons
   );

/** creates the master var and initializes the vardata */
extern
SCIP_RETCODE GCGcreateMasterVar(
   SCIP* scip,
   SCIP* pricingscip,
   SCIP_VAR** newvar,
   char* varname,
   SCIP_Real objcoeff,
   SCIP_VARTYPE vartype,
   SCIP_Bool solisray,
   int prob,
   int nsolvars,
   SCIP_Real* solvals,
   SCIP_VAR** solvars
   );

/** creates initial master variables and the vardata */
SCIP_RETCODE GCGcreateInitialMasterVar(
   SCIP* scip,
   SCIP_VAR* var,
   SCIP_VAR** newvar
   );


#ifdef __cplusplus
}
#endif

#endif
