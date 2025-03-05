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

/**@file   gcgcol.c
 * @brief  methods for working with gcg column structure
 * @author Jonas Witt
 *
 * Various methods to work with gcg column structure
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/pub_gcgcol.h"

#include "gcg/gcg.h"
#include "scip/def.h"
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "gcg/pricer_gcg.h"
#include "gcg/sepa_original.h"
#include "gcg/struct_gcgcol.h"

#include <assert.h>

/** create a gcg column */
SCIP_RETCODE GCGcreateGcgCol(
   SCIP*                pricingprob,        /**< SCIP data structure (pricing problem) */
   GCG_COL**            gcgcol,             /**< pointer to store gcg column */
   int                  probnr,             /**< number of corresponding pricing problem */
   SCIP_VAR**           vars,               /**< (sorted) array of variables of corresponding pricing problem */
   SCIP_Real*           vals,               /**< array of solution values (belonging to vars) */
   int                  nvars,              /**< number of variables */
   SCIP_Bool            isray,              /**< is the column a ray? */
   SCIP_Real            redcost             /**< last known reduced cost */
   )
{
   int i;
   int nnonz;
   int size;

   /* WARNING: this function has to be threadsafe!*/

   size = SCIPcalcMemGrowSize(pricingprob, nvars);
   SCIP_CALL( SCIPallocBlockMemory(pricingprob, gcgcol) );
   SCIP_CALL( SCIPallocBlockMemoryArray(pricingprob, &((*gcgcol)->vars), size) );
   SCIP_CALL( SCIPallocBlockMemoryArray(pricingprob, &((*gcgcol)->vals), size) );

   (*gcgcol)->pricingprob = pricingprob;
   (*gcgcol)->probnr = probnr;
   (*gcgcol)->isray = isray;
   (*gcgcol)->redcost = redcost;
   (*gcgcol)->age = 0;
   (*gcgcol)->mastercoefs = NULL;
   (*gcgcol)->originalsepamastercuts = NULL;
   (*gcgcol)->linkvars = NULL;
   (*gcgcol)->nmastercoefs = 0;
   (*gcgcol)->noriginalsepamastercuts = 0;
   (*gcgcol)->maxmastercoefs = 0;
   (*gcgcol)->maxoriginalsepamastercuts = 0;
   (*gcgcol)->extendedmasterconscoefs = NULL;
   (*gcgcol)->extendedmasterconsbounds = NULL;
   (*gcgcol)->nextendedmasterconss = 0;
   (*gcgcol)->nlinkvars = 0;
   (*gcgcol)->initcoefs = FALSE;
   (*gcgcol)->pos = -1;
   (*gcgcol)->maxvars = size;

   nnonz = 0;
   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* origvar;
      SCIP_Real scalar;
      SCIP_Real constant;
      SCIP_Real origval;

      scalar = 1.0;
      constant = 0.0;

      origvar = vars[i];

      /* todo: capture vars? */
      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

      assert( !SCIPisZero(pricingprob, scalar) );

      origval = (vals[i] - constant) / scalar;

      /* round origval if possible to avoid numerical troubles */
      if( SCIPvarIsIntegral(origvar) && SCIPisFeasIntegral(pricingprob, origval) )
         origval = SCIPround(pricingprob, origval);

      if( SCIPisZero(pricingprob, origval) )
         continue;

      assert((GCGvarIsPricing(origvar) && GCGpricingVarGetNOrigvars(origvar) > 0 && GCGpricingVarGetOrigvars(origvar)[0] != NULL) || GCGvarIsInferredPricing(origvar));

      SCIPcaptureVar(pricingprob, origvar);
      (*gcgcol)->vars[nnonz] = origvar;
      (*gcgcol)->vals[nnonz] = origval;
      ++nnonz;
   }

   (*gcgcol)->nvars = nnonz;

   /* sort vars and vals array w.r.t. variable index */
   SCIPsortPtrReal((void**)(*gcgcol)->vars, (double*)(*gcgcol)->vals, SCIPvarComp, nnonz);

#ifndef NDEBUG
   for( i = 1 ; i < (*gcgcol)->nvars; ++i )
   {
      assert( SCIPvarCompare((*gcgcol)->vars[i-1], (*gcgcol)->vars[i]) < 0 );
   }
#endif
   return SCIP_OKAY;
}

/** free a gcg column */
void GCGfreeGcgCol(
   GCG_COL**            gcgcol              /**< pointer to store gcg column */
   )
{
   int i;
   assert(gcgcol != NULL);
   assert(*gcgcol != NULL);

   /* WARNING: this function has to be threadsafe!*/

   for( i = 0; i < (*gcgcol)->nvars; ++i )
      SCIPreleaseVar((*gcgcol)->pricingprob, &(*gcgcol)->vars[i]);
   assert((*gcgcol)->nvars == 0 || (*gcgcol)->vars != NULL);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->vars, (*gcgcol)->maxvars);
   assert((*gcgcol)->nvars == 0 || (*gcgcol)->vals != NULL);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->vals, (*gcgcol)->maxvars);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->mastercoefs, (*gcgcol)->maxmastercoefs);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->linkvars, (*gcgcol)->maxlinkvars);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->originalsepamastercuts, (*gcgcol)->maxoriginalsepamastercuts);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->extendedmasterconscoefs, (*gcgcol)->nextendedmasterconss);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->extendedmasterconsbounds, (*gcgcol)->nextendedmasterconss);
   SCIPfreeBlockMemory((*gcgcol)->pricingprob, gcgcol);
}

/** create a gcg column from a solution to a pricing problem */
SCIP_RETCODE GCGcreateGcgColFromSol(
   SCIP*                pricingprob,        /**< SCIP data structure (pricing problem) */
   SCIP*                subproblem,         /**< SCIP data structure that contains the actual solution (if NULL pricingprob will be used) */
   SCIP_HASHMAP*        varmap,             /**< mapping of pricingprob vars to subproblem vars (can be NULL if subproblem is NULL) */
   GCG_COL**            gcgcol,             /**< pointer to store gcg column */
   int                  prob,               /**< number of corresponding pricing problem */
   SCIP_SOL*            sol,                /**< solution of pricing problem with index prob */
   SCIP_Bool            isray,              /**< is column a ray? */
   SCIP_Real            redcost             /**< last known reduced cost */
   )
{
   SCIP* solprob;
   SCIP_VAR** solvars;
   SCIP_VAR** colvars;

   SCIP_Real* colvals;

   int nsolvars;
   int ncolvars;

   int i;

   if ( subproblem == NULL )
      solprob = pricingprob;
   else
      solprob = subproblem;

   solvars = SCIPgetOrigVars(pricingprob);
   nsolvars = SCIPgetNOrigVars(pricingprob);

   SCIP_CALL( SCIPallocBufferArray(pricingprob, &colvars, nsolvars) );
   SCIP_CALL( SCIPallocBufferArray(pricingprob, &colvals, nsolvars) );

   ncolvars = 0;

   for( i = 0; i < nsolvars; ++i )
   {
      SCIP_VAR* solvar;
      SCIP_Real solval;

      solvar = solvars[i];
      if ( varmap != NULL )
         solval = SCIPgetSolVal(subproblem, sol, SCIPhashmapGetImage(varmap, solvar));
      else
         solval = SCIPgetSolVal(pricingprob, sol, solvar);

      /* round solval if possible to avoid numerical troubles */
      if( SCIPvarIsIntegral(solvar) && SCIPisFeasIntegral(solprob, solval) )
         solval = SCIPround(solprob, solval);

      if( SCIPisZero(solprob, solval) )
      {
         continue;
      }

      colvars[ncolvars] = solvar;
      colvals[ncolvars] = solval;
      ++ncolvars;
   }

   SCIP_CALL( GCGcreateGcgCol(pricingprob, gcgcol, prob, colvars, colvals, ncolvars, isray, redcost) );

   SCIPfreeBufferArray(pricingprob, &colvals);
   SCIPfreeBufferArray(pricingprob, &colvars);

   return SCIP_OKAY;
}

/** comparison method for sorting gcg columns by non-decreasing reduced cost */
SCIP_DECL_SORTPTRCOMP(GCGcolCompRedcost)
{
   SCIP_Real redcost1;
   SCIP_Real redcost2;

   redcost1 = GCGcolGetRedcost((GCG_COL*) elem1);
   redcost2 = GCGcolGetRedcost((GCG_COL*) elem2);

   if( redcost1 < redcost2 )
      return -1;
   else if( redcost1 > redcost2 )
      return +1;
   else
      return 0;
}

/** comparison method for sorting gcg columns by non-increasing age */
SCIP_DECL_SORTPTRCOMP(GCGcolCompAge)
{
   int age1;
   int age2;

   age1 = GCGcolGetAge((GCG_COL*) elem1);
   age2 = GCGcolGetAge((GCG_COL*) elem2);

   if( age1 < age2 )
      return +1;
   else if( age1 > age2 )
      return -1;
   else
      return 0;
}

/** comparison method for gcg columns. Returns TRUE iff columns are equal */
SCIP_Bool GCGcolIsEq(
   GCG_COL*             gcgcol1,
   GCG_COL*             gcgcol2
   )
{
   SCIP* pricingprob;

   SCIP_VAR** vars1;
   SCIP_VAR** vars2;

   SCIP_Real* vals1;
   SCIP_Real* vals2;

   int nvars1;
   int nvars2;
   int probnr1;
   int probnr2;

   int i;


   probnr1 = GCGcolGetProbNr(gcgcol1);
   probnr2 = GCGcolGetProbNr(gcgcol2);

   if( probnr1 != probnr2 )
      return FALSE;

   nvars1 = GCGcolGetNVars(gcgcol1);
   nvars2 = GCGcolGetNVars(gcgcol2);

   if( nvars1 != nvars2 )
      return FALSE;

   pricingprob = GCGcolGetPricingProb(gcgcol1);
   vars1 = GCGcolGetVars(gcgcol1);
   vars2 = GCGcolGetVars(gcgcol2);

   vals1 = GCGcolGetVals(gcgcol1);
   vals2 = GCGcolGetVals(gcgcol2);

   for( i = 0; i < nvars1; ++i )
   {
      if( vars1[i] != vars2[i] || !SCIPisEQ(pricingprob, vals1[i], vals2[i]) )
      {
         return FALSE;
      }
   }

   return TRUE;

}

/** get pricing problem index of gcg column */
int GCGcolGetProbNr(
   GCG_COL*             gcgcol
   )
{
   return gcgcol->probnr;
}

/** get pricing problem of gcg column */
SCIP* GCGcolGetPricingProb(
   GCG_COL*             gcgcol
   )
{
   return gcgcol->pricingprob;
}

/** get variables of gcg column */
SCIP_VAR** GCGcolGetVars(
   GCG_COL*             gcgcol
   )
{
   return gcgcol->vars;
}

/** get values of gcg column */
SCIP_Real* GCGcolGetVals(
   GCG_COL*             gcgcol
   )
{
   return gcgcol->vals;
}

/** get number of variables of gcg column */
int GCGcolGetNVars(
   GCG_COL*             gcgcol
   )
{
   return gcgcol->nvars;
}

/** is gcg column a ray? */
SCIP_Bool GCGcolIsRay(
   GCG_COL*             gcgcol
   )
{
   return gcgcol->isray;
}

/** get reduced cost of gcg column */
SCIP_Real GCGcolGetRedcost(
   GCG_COL*             gcgcol
   )
{
   return gcgcol->redcost;
}

/** get age of gcg column */
int GCGcolGetAge(
   GCG_COL*             gcgcol
   )
{
   return gcgcol->age;
}

/** update reduced cost of variable and increase age */
void GCGcolUpdateRedcost(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   SCIP_Real            redcost,            /**< new reduced cost */
   SCIP_Bool            growold             /**< change age depending on reduced cost? */
   )
{
   gcgcol->redcost = redcost;

   if( !growold )
      return;

   if( !SCIPisNegative(gcgcol->pricingprob, redcost) )
      ++(gcgcol->age);
   else
      gcgcol->age = 0;
}

/** get master coefficients of column */
SCIP_Real* GCGcolGetMastercoefs(
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   return gcgcol->mastercoefs;
}

/** get number of master coefficients of column */
int GCGcolGetNMastercoefs(
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   return gcgcol->nmastercoefs;
}

/** set master coefficients information of column */
SCIP_RETCODE GCGcolSetMastercoefs(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   SCIP_Real*           mastercoefs,        /**< array of master coefficients */
   int                  nmastercoefs        /**< new number of master coefficients */
   )
{
   int i;

   SCIPdebugMessage("Col set master coefs\n");
   assert(gcgcol->nmastercoefs == 0);
   if( nmastercoefs == 0 )
      return SCIP_OKAY;

   gcgcol->maxmastercoefs = SCIPcalcMemGrowSize(gcgcol->pricingprob, nmastercoefs);
   SCIP_CALL( SCIPallocBlockMemoryArray(gcgcol->pricingprob, &(gcgcol->mastercoefs), gcgcol->maxmastercoefs) );

   for( i = 0; i < nmastercoefs; ++i )
   {
      SCIP_Real coef = mastercoefs[i];
      gcgcol->mastercoefs[i] = coef;
   }

   gcgcol->nmastercoefs = nmastercoefs;

   return SCIP_OKAY;
}

/** set norm of column */
void GCGcolSetNorm(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   SCIP_Real            norm                /**< norm of column */
   )
{
   gcgcol->norm = norm;
}

/** get norm of column */
void GCGcolComputeNorm(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   int i;
   SCIP_Real norm;

   assert(scip != NULL);
   assert(gcgcol != NULL);

   norm = 0.0;
   /* compute scalar of master values of gcg columns */
   for( i = 0; i < gcgcol->nmastercoefs; ++i )
   {
      if( !SCIPisZero(scip, gcgcol->mastercoefs[i]) )
         norm += SQR(gcgcol->mastercoefs[i]);
   }

   for( i = 0; i < gcgcol->noriginalsepamastercuts; ++i )
   {
      if( !SCIPisZero(scip, gcgcol->originalsepamastercuts[i]) )
         norm += SQR(gcgcol->originalsepamastercuts[i]);
   }

   for( i = 0; i < gcgcol->nextendedmasterconss; ++i )
   {
      if( !SCIPisZero(scip, gcgcol->extendedmasterconscoefs[i]) )
         norm += SQR(gcgcol->extendedmasterconscoefs[i]);
   }

   for( i = 0; i < gcgcol->nlinkvars; ++i )
   {
      if( !SCIPisZero(scip, gcgcol->vals[gcgcol->linkvars[i]]) )
         norm += gcgcol->vals[gcgcol->linkvars[i]];
   }

   /* consider convexity constraint */
   norm += 1.0;

   gcgcol->norm = norm;
}

/** set master coefficients of column as initialized */
SCIP_RETCODE GCGcolSetInitializedCoefs(
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   assert(!gcgcol->initcoefs);
   gcgcol->initcoefs = TRUE;
   return SCIP_OKAY;
}

/** return if master coefficients of column have been initialized */
SCIP_Bool GCGcolGetInitializedCoefs(
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   return gcgcol->initcoefs;
}

/** get master coefficients of column */
int* GCGcolGetLinkvars(
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   return gcgcol->linkvars;
}

/** get number of master coefficients of column */
int GCGcolGetNLinkvars(
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   return gcgcol->nlinkvars;
}

/** set master coefficients information of column */
SCIP_RETCODE GCGcolSetLinkvars(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   int*                 linkvars,           /**< array of linking variable indices for gcgcol->var */
   int                  nlinkvars           /**< number of linking variables in gcgcol->var */
   )
{
   int i;

   assert(gcgcol->nlinkvars == 0);

   gcgcol->maxlinkvars = SCIPcalcMemGrowSize(gcgcol->pricingprob, nlinkvars);
   SCIP_CALL( SCIPallocBlockMemoryArray(gcgcol->pricingprob, &(gcgcol->linkvars), gcgcol->maxlinkvars) );

   for( i = 0; i < nlinkvars; ++i )
   {
      gcgcol->linkvars[i] = linkvars[i];
   }

   gcgcol->nlinkvars = nlinkvars;

   return SCIP_OKAY;
}

/** get original separator cut coefficients of column in master problem */
SCIP_Real* GCGcolGetOriginalSepaMastercuts(
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   return gcgcol->originalsepamastercuts;
}

/** get number of original separator cut coefficients of column in master problem */
int GCGcolGetNOriginalSepaMastercuts(
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   return gcgcol->noriginalsepamastercuts;
}

/** get extended master cons coefficients of column in the master problem */
GCG_EXPORT
SCIP_Real* GCGcolGetExtendedmastercons(
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   return gcgcol->extendedmasterconscoefs;
}

/** get extended master cons bounds */
static
SCIP_Real* colGetExtendedmasterconsBounds(
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   return gcgcol->extendedmasterconsbounds;
}

/** get number of extended master cons coefficients of column in the master problem */
GCG_EXPORT
int GCGcolGetNExtendedmasterconss(
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   return gcgcol->nextendedmasterconss;
}

/** get norm of column */
SCIP_Real GCGcolGetNorm(
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   return gcgcol->norm;
}

/** update original separator cut coefficients information of column in the master problem */
SCIP_RETCODE GCGcolUpdateOriginalSepaMastercuts(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   SCIP_Real*           neworiginalsepamastercuts,/**< pointer to new array of extended master cons coefficients */
   int                  nneworiginalsepamastercuts/**< new number of extended master cons coefficients */
   )
{
   int i;
   int newsize;

   i = gcgcol->noriginalsepamastercuts + nneworiginalsepamastercuts;
   newsize = SCIPcalcMemGrowSize(gcgcol->pricingprob, i);
   if( i > gcgcol->maxoriginalsepamastercuts )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(GCGcolGetPricingProb(gcgcol), &(gcgcol->originalsepamastercuts),
            gcgcol->maxoriginalsepamastercuts, newsize) );
   }

   gcgcol->maxoriginalsepamastercuts = newsize;

   for( i = 0; i < nneworiginalsepamastercuts; ++i )
   {
      gcgcol->originalsepamastercuts[gcgcol->noriginalsepamastercuts] = neworiginalsepamastercuts[i];
      ++(gcgcol->noriginalsepamastercuts);
   }

   return SCIP_OKAY;
}

/** set extended master cons coefficients information of column in the master problem
 * @note the arrays will be freed by the column, they must be allocated using the pricingscip the column belongs to
 */
GCG_EXPORT
SCIP_RETCODE GCGcolSetExtendedmasterconss(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   SCIP_Real*           extendedmasterconss,  /**< pointer to array of extended master cons coefficients */
   SCIP_Real*           extendedmasterconsbounds,/**< pointer to array of extended master cons bounds */
   int                  nextendedmasterconss  /**< number of extended master cons coefficients */
   )
{
   SCIPdebugMessage("Col set master coefs\n");

   if( gcgcol->nextendedmasterconss > 0 )
   {
      SCIPfreeBlockMemoryArrayNull(gcgcol->pricingprob, &(gcgcol->extendedmasterconscoefs), gcgcol->nextendedmasterconss);
      SCIPfreeBlockMemoryArrayNull(gcgcol->pricingprob, &(gcgcol->extendedmasterconsbounds), gcgcol->nextendedmasterconss);
   }

   gcgcol->nextendedmasterconss = nextendedmasterconss;

   assert(gcgcol->extendedmasterconscoefs == NULL);
   assert(gcgcol->extendedmasterconsbounds == NULL);

   if( nextendedmasterconss == 0 )
      return SCIP_OKAY;

   gcgcol->extendedmasterconscoefs = extendedmasterconss;
   gcgcol->extendedmasterconsbounds = extendedmasterconsbounds;

   return SCIP_OKAY;
}

/** return solution value of variable in gcg column */
SCIP_Real GCGcolGetSolVal(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_COL*             gcgcol,             /**< gcg column */
   SCIP_VAR*            var                 /**< variable */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int nvars;
   int pos;
   SCIP_Bool found;

   vars = gcgcol->vars;
   vals = gcgcol->vals;
   nvars = gcgcol->nvars;

   found = SCIPsortedvecFindPtr((void**) vars, SCIPvarComp, (void*) var, nvars, &pos);

   if( !found )
   {
      return 0.0;
   }

   return vals[pos];
}

/** returns true if the gcg column knows the solution value of the variable */
GCG_EXPORT
SCIP_Bool GCGcolKnowsSolVar(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_COL*             gcgcol,             /**< gcg column */
   SCIP_VAR*            var                 /**< variable */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int pos;

   vars = gcgcol->vars;
   nvars = gcgcol->nvars;

   return SCIPsortedvecFindPtr((void**) vars, SCIPvarComp, (void*) var, nvars, &pos);
}

/** returns whether the col's age exceeds the age limit */
SCIP_Bool GCGcolIsAged(
   GCG_COL*             col,                /**< col to check */
   int                   agelimit            /**< maximum age a col can reach before it is deleted from the pool, or -1 */
   )
{
   assert(col != NULL);

   return (agelimit >= 0 && col->age > agelimit);
}

/** compute parallelism of column to dual objective */
SCIP_Real GCGcolComputeDualObjPara(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_COL*             gcgcol              /**< gcg column */
   )
{
   SCIP_Real para;
   int i;
   SCIP_CONS** masterconss;
   SCIP_ROW** originalsepamastercuts;
   SCIP_Real* extendedmasterconsbounds;
   SCIP_Real dualobjnorm;

   assert(scip != NULL);
   assert(gcgcol != NULL);

   masterconss = GCGgetMasterConss(GCGmasterGetOrigprob(scip));
   originalsepamastercuts = GCGsepaGetOriginalSepaMastercuts(scip);
   extendedmasterconsbounds = colGetExtendedmasterconsBounds(gcgcol);

   para = 0.0;

   dualobjnorm = 0.0;

   /* compute scalar of master values of gcg columns */
   for( i = 0; i < gcgcol->nmastercoefs; ++i )
   {
      SCIP_Real lhs;
      SCIP_Real rhs;

      lhs = SCIPgetLhsLinear(scip, masterconss[i]);
      rhs = SCIPgetRhsLinear(scip, masterconss[i]);

      if( !SCIPisInfinity(scip, -lhs) )
      {
         dualobjnorm += SQR(lhs);

         if( SCIPisPositive(scip, gcgcol->mastercoefs[i]) )
            para += gcgcol->mastercoefs[i] * lhs;
      }
      else if( !SCIPisInfinity(scip, rhs) )
      {
         dualobjnorm += SQR(rhs);

         if(SCIPisNegative(scip, gcgcol->mastercoefs[i] ) )
            para += gcgcol->mastercoefs[i] * rhs;
      }
   }

   for( i = 0; i < gcgcol->noriginalsepamastercuts; ++i )
   {
      SCIP_Real lhs;
      SCIP_Real rhs;

      if( !SCIProwIsInLP(originalsepamastercuts[i]) )
         continue;

      lhs = SCIProwGetLhs(originalsepamastercuts[i]);
      rhs = SCIProwGetRhs(originalsepamastercuts[i]);

      if( !SCIPisInfinity(scip, -lhs))
      {
         dualobjnorm += SQR(lhs);

         if( SCIPisPositive(scip, gcgcol->originalsepamastercuts[i]) )
            para += gcgcol->originalsepamastercuts[i] * lhs;
      }
      else if( !SCIPisInfinity(scip, rhs) )
      {
         dualobjnorm += SQR(rhs);

         if(SCIPisNegative(scip, gcgcol->originalsepamastercuts[i] ) )
            para += gcgcol->originalsepamastercuts[i] * rhs;
      }
   }

   for( i = 0; i < gcgcol->nextendedmasterconss; ++i )
   {
      SCIP_Real bound = extendedmasterconsbounds[i];

      if( SCIPisInfinity(scip, ABS(bound)) )
      {
         return SCIP_ERROR;
      }

      dualobjnorm += SQR(bound);

      if( SCIPisPositive(scip, gcgcol->extendedmasterconscoefs[i]) )
         para += gcgcol->extendedmasterconscoefs[i] * bound;
   }

   for( i = 0; i < GCGgetNPricingprobs(GCGmasterGetOrigprob(scip)); ++i )
      dualobjnorm += SQR(GCGgetNIdenticalBlocks(GCGmasterGetOrigprob(scip), i));

   para += SQR(GCGgetNIdenticalBlocks(GCGmasterGetOrigprob(scip), gcgcol->probnr));

   assert(!SCIPisInfinity(scip, ABS(para)));

   dualobjnorm = sqrt(dualobjnorm);
   assert(!SCIPisInfinity(scip, dualobjnorm));
   assert(SCIPisPositive(scip, dualobjnorm));
   assert(SCIPisPositive(scip, gcgcol->norm));

   para = para / (dualobjnorm * gcgcol->norm);

   return para;
}

/** compute orthogonality of two gcg columns */
SCIP_Real GCGcolComputeOrth(
   SCIP*                scip,               /**< SCIP data structure */
   GCG_COL*             gcgcol1,            /**< first gcg column */
   GCG_COL*             gcgcol2             /**< second gcg column */
   )
{
   int i;
   int j;
   SCIP_Real para = 0.0;
   SCIP_Real norm1 = 0.0;
   SCIP_Real norm2 = 0.0;

   assert(scip != NULL);
   assert(gcgcol1 != NULL);
   assert(gcgcol2 != NULL);

   /* compute scalar of master values of gcg columns */
   for( i = 0; i < gcgcol1->nmastercoefs; ++i )
   {
      if( SCIPisPositive(scip, gcgcol1->mastercoefs[i] * gcgcol2->mastercoefs[i]) )
         para += gcgcol1->mastercoefs[i] * gcgcol2->mastercoefs[i];

      if( SCIPisPositive(scip, gcgcol1->mastercoefs[i]) )
         norm1 += SQR(gcgcol1->mastercoefs[i]);
      if( SCIPisPositive(scip, gcgcol2->mastercoefs[i]) )
         norm2 += SQR(gcgcol2->mastercoefs[i]);
   }

   for( i = 0; i < gcgcol1->noriginalsepamastercuts; ++i )
   {
      if( SCIPisPositive(scip, gcgcol1->originalsepamastercuts[i] * gcgcol2->originalsepamastercuts[i]) )
         para += gcgcol1->originalsepamastercuts[i] * gcgcol2->originalsepamastercuts[i];

      if( SCIPisPositive(scip, gcgcol1->originalsepamastercuts[i]) )
         norm1 += SQR(gcgcol1->originalsepamastercuts[i]);
      if( SCIPisPositive(scip, gcgcol2->originalsepamastercuts[i]) )
         norm2 += SQR(gcgcol2->originalsepamastercuts[i]);
   }

   for( i = 0; i < gcgcol1->nextendedmasterconss; ++i )
   {
      if( SCIPisPositive(scip, gcgcol1->extendedmasterconscoefs[i] * gcgcol2->extendedmasterconscoefs[i]) )
         para += gcgcol1->extendedmasterconscoefs[i] * gcgcol2->extendedmasterconscoefs[i];

      if( SCIPisPositive(scip, gcgcol1->extendedmasterconscoefs[i]) )
         norm1 += SQR(gcgcol1->extendedmasterconscoefs[i]);
      if( SCIPisPositive(scip, gcgcol2->extendedmasterconscoefs[i]) )
         norm2 += SQR(gcgcol2->extendedmasterconscoefs[i]);
   }

   for( i = 0; i < gcgcol1->nlinkvars; ++i )
   {
      SCIP_VAR* linkvar1;
      SCIP_Real linkval1;
      linkvar1 = gcgcol1->vars[gcgcol1->linkvars[i]];
      linkval1 = gcgcol1->vals[gcgcol1->linkvars[i]];

      norm1 += SQR(linkval1);

      for( j = 0; j < gcgcol2->nlinkvars; ++j )
      {
         SCIP_VAR* linkvar2;
         SCIP_Real linkval2;
         linkvar2 = gcgcol2->vars[gcgcol2->linkvars[j]];
         linkval2 = gcgcol2->vals[gcgcol2->linkvars[j]];

         if( linkvar1 == linkvar2 )
         {
            para += linkval1 * linkval2;
            break;
         }
      }
   }

   for( i = 0; i < gcgcol2->nlinkvars; ++i )
   {
      SCIP_Real linkval2;

      linkval2 = gcgcol2->vals[gcgcol2->linkvars[i]];

      norm2 += SQR(linkval2);
   }


   /* scalar for convexitiy constraints */
   if( gcgcol1->probnr == gcgcol2->probnr )
      para *= 1.0;

   norm1 *= 1.0;
   norm2 *= 1.0;

   norm1 = sqrt(norm1);
   norm2 = sqrt(norm2);

   assert(SCIPisPositive(scip, norm1) && SCIPisPositive(scip, norm2));

   para = para/(norm1*norm2);

   return 1.0 - para;
}

SCIP_DECL_HASHGETKEY(GCGhashGetKeyCol)
{  /*lint --e{715}*/
   assert(elem != NULL);

   /* the key of a col is the col itself */
   return elem;
}

SCIP_DECL_HASHKEYEQ(GCGhashKeyEqCol)
{  /*lint --e{715}*/
   return GCGcolIsEq((GCG_COL*)key1, (GCG_COL*)key2);
}

SCIP_DECL_HASHKEYVAL(GCGhashKeyValCol)
{  /*lint --e{715}*/
   GCG_COL* col;
   unsigned int keyval;
   int minindex;
   int maxindex;

   col = (GCG_COL*)key;
   assert(col != NULL);

   /* TODO: this hash function does not respect tolerances (except the hard coded of SCIPrealHashCode)
    * but it seems that SCIP does it the same way */
   if( col->nvars > 0 )
   {
      minindex = SCIPvarGetIndex(col->vars[0]);
      maxindex = SCIPvarGetIndex(col->vars[col->nvars-1]);
   }
   else
   {
      minindex = INT_MAX;
      maxindex = INT_MAX;
   }
   assert(minindex <= maxindex);
   keyval = SCIPhashSeven(col->probnr, col->nvars, col->isray,
         SCIPrealHashCode(col->nvars > 0 ? col->vals[0] : 0.0), minindex,
         SCIPrealHashCode(col->nvars > 0 ? col->vals[col->nvars-1] : 0.0), maxindex);

   return keyval;
}
