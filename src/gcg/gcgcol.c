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
#include "gcg/struct_mastersepacut.h"

#include <assert.h>

/** create a gcg column */
SCIP_RETCODE GCGcreateGcgCol(
   GCG*                 gcg,                /**< GCG data structure */
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
   int ninferredpricingvars;

   /* WARNING: this function has to be threadsafe!*/

   ninferredpricingvars = GCGcountInferredCoefPricingVars(vars, nvars);

   SCIP_CALL( SCIPallocBlockMemory(pricingprob, gcgcol) );
   (*gcgcol)->maxvars = SCIPcalcMemGrowSize(pricingprob, nvars - ninferredpricingvars);
   SCIP_CALL( SCIPallocBlockMemoryArray(pricingprob, &((*gcgcol)->vars), (*gcgcol)->maxvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(pricingprob, &((*gcgcol)->vals), (*gcgcol)->maxvars) );

   if( ninferredpricingvars > 0 )
   {
      (*gcgcol)->maxinferredpricingvars = SCIPcalcMemGrowSize(pricingprob, ninferredpricingvars);
      SCIP_CALL( SCIPallocBlockMemoryArray(pricingprob, &((*gcgcol)->inferredpricingvars), (*gcgcol)->maxinferredpricingvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(pricingprob, &((*gcgcol)->inferredpricingvals), (*gcgcol)->maxinferredpricingvars) );
   }
   else
   {
      (*gcgcol)->inferredpricingvars = NULL;
      (*gcgcol)->inferredpricingvals = NULL;
      (*gcgcol)->maxinferredpricingvars = 0;
   }

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
   (*gcgcol)->nlinkvars = 0;
   (*gcgcol)->initcoefs = FALSE;

   /* data structures to store the cut coefficients for active master separator cuts */
   (*gcgcol)->sepamastercutcoeffs = NULL;
   (*gcgcol)->nsepamastercutcoeffs = 0;
   (*gcgcol)->sepamastercutscoeffssize = 0;

   (*gcgcol)->pos = -1;


   (*gcgcol)->nvars = 0;
   (*gcgcol)->ninferredpricingvars = 0;
   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* origvar;
      SCIP_Real scalar;
      SCIP_Real constant;
      SCIP_Real origval;

      scalar = 1.0;
      constant = 0.0;

      origvar = vars[i];

      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

      assert( !SCIPisZero(pricingprob, scalar) );

      origval = (vals[i] - constant) / scalar;

      /* round origval if possible to avoid numerical troubles */
      if( SCIPvarIsIntegral(origvar) && SCIPisFeasIntegral(pricingprob, origval) )
         origval = SCIPround(pricingprob, origval);

      if( SCIPisZero(pricingprob, origval) )
         continue;

      if( GCGvarIsInferredPricing(origvar) && GCGinferredPricingVarIsCoefVar(origvar) )
      {
         assert((*gcgcol)->inferredpricingvars != NULL);
         assert((*gcgcol)->ninferredpricingvars < (*gcgcol)->maxinferredpricingvars);
         (*gcgcol)->inferredpricingvars[(*gcgcol)->ninferredpricingvars] = origvar;
         (*gcgcol)->inferredpricingvals[(*gcgcol)->ninferredpricingvars] = origval;
         ++(*gcgcol)->ninferredpricingvars;
         SCIPcaptureVar(pricingprob, origvar);
      }
      else if( GCGvarIsPricing(origvar) )
      {
         assert(GCGvarIsPricing(origvar) && GCGpricingVarGetNOrigvars(origvar) > 0 && GCGpricingVarGetOrigvars(origvar)[0] != NULL);
         assert((*gcgcol)->nvars < (*gcgcol)->maxvars);
         (*gcgcol)->vars[(*gcgcol)->nvars] = origvar;
         (*gcgcol)->vals[(*gcgcol)->nvars] = origval;
         ++(*gcgcol)->nvars;
         SCIPcaptureVar(pricingprob, origvar);
      }
   }
   assert(ninferredpricingvars == (*gcgcol)->ninferredpricingvars);

   /* sort vars and vals array w.r.t. variable index */
   SCIPsortPtrReal((void**)(*gcgcol)->vars, (double*)(*gcgcol)->vals, SCIPvarComp, (*gcgcol)->nvars);
   if( (*gcgcol)->inferredpricingvars != NULL )
   {
      SCIPsortPtrReal((void**)(*gcgcol)->inferredpricingvars, (double*)(*gcgcol)->inferredpricingvals,
         SCIPvarComp, (*gcgcol)->ninferredpricingvars);
   }

#ifndef NDEBUG
   for( i = 1 ; i < (*gcgcol)->nvars; ++i )
   {
      assert( SCIPvarCompare((*gcgcol)->vars[i-1], (*gcgcol)->vars[i]) < 0 );
   }
   for( i = 1 ; i < (*gcgcol)->ninferredpricingvars; ++i )
   {
      assert( SCIPvarCompare((*gcgcol)->inferredpricingvars[i-1], (*gcgcol)->inferredpricingvars[i]) < 0 );
   }
#endif
   return SCIP_OKAY;
}

/** free a gcg column */
SCIP_RETCODE GCGfreeGcgCol(
   GCG_COL**            gcgcol              /**< pointer to store gcg column */
   )
{
   int i;
   assert(gcgcol != NULL);
   assert(*gcgcol != NULL);

   /* WARNING: this function has to be threadsafe!*/

   for( i = 0; i < (*gcgcol)->nvars; ++i )
      SCIP_CALL( SCIPreleaseVar((*gcgcol)->pricingprob, &(*gcgcol)->vars[i]) );
   assert((*gcgcol)->nvars == 0 || (*gcgcol)->vars != NULL);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->vars, (*gcgcol)->maxvars);
   assert((*gcgcol)->nvars == 0 || (*gcgcol)->vals != NULL);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->vals, (*gcgcol)->maxvars);
   for( i = 0; i < (*gcgcol)->ninferredpricingvars; ++i )
      SCIP_CALL( SCIPreleaseVar((*gcgcol)->pricingprob, &(*gcgcol)->inferredpricingvars[i]) );
   assert((*gcgcol)->ninferredpricingvars == 0 || (*gcgcol)->inferredpricingvars != NULL);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->inferredpricingvars, (*gcgcol)->maxinferredpricingvars);
   assert((*gcgcol)->ninferredpricingvars == 0 || (*gcgcol)->inferredpricingvals != NULL);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->inferredpricingvals, (*gcgcol)->maxinferredpricingvars);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->mastercoefs, (*gcgcol)->maxmastercoefs);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->linkvars, (*gcgcol)->maxlinkvars);
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &(*gcgcol)->originalsepamastercuts, (*gcgcol)->maxoriginalsepamastercuts);

   /* free data structures for storing coefficients of mastercuts generated by separators */
   SCIPfreeBlockMemoryArrayNull((*gcgcol)->pricingprob, &((*gcgcol)->sepamastercutcoeffs), (*gcgcol)->sepamastercutscoeffssize);
   (*gcgcol)->sepamastercutscoeffssize =0;
   (*gcgcol)->nsepamastercutcoeffs = 0;

   SCIPfreeBlockMemory((*gcgcol)->pricingprob, gcgcol);
   return SCIP_OKAY;
}

/** create a gcg column from a solution to a pricing problem */
SCIP_RETCODE GCGcreateGcgColFromSol(
   GCG*                 gcg,                /**< GCG data structure */
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

   SCIP_CALL( GCGcreateGcgCol(gcg, pricingprob, gcgcol, prob, colvars, colvals, ncolvars, isray, redcost) );

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

   int i;

   if( gcgcol1->probnr != gcgcol2->probnr )
      return FALSE;

   if( gcgcol1->nvars != gcgcol2->nvars )
      return FALSE;
   
   if( gcgcol1->ninferredpricingvars != gcgcol2->ninferredpricingvars )
      return FALSE;

   pricingprob = GCGcolGetPricingProb(gcgcol1);

   vars1 = GCGcolGetVars(gcgcol1);
   vars2 = GCGcolGetVars(gcgcol2);

   vals1 = GCGcolGetVals(gcgcol1);
   vals2 = GCGcolGetVals(gcgcol2);

   for( i = 0; i < gcgcol1->nvars; ++i )
   {
      if( vars1[i] != vars2[i] || !SCIPisEQ(pricingprob, vals1[i], vals2[i]) )
      {
         return FALSE;
      }
   }

   vars1 = gcgcol1->inferredpricingvars;
   vars2 = gcgcol2->inferredpricingvars;

   vals1 = gcgcol1->inferredpricingvals;
   vals2 = gcgcol2->inferredpricingvals;

   for( i = 0; i < gcgcol1->ninferredpricingvars; ++i )
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
SCIP_RETCODE GCGcolComputeNorm(
   GCG*                 gcg,                /**< GCG data structure */
   GCG_COL*             gcgcol              /**< gcg column structure */
   )
{
   SCIP* scip;
   int i;
   int j;
   SCIP_Real norm;

   assert(gcg != NULL);
   assert(gcgcol != NULL);

   scip = GCGgetMasterprob(gcg);
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

   /* consider coefficients for cuts generated by separators */
   for( j = 0; j < gcgcol->nsepamastercutcoeffs; j++ )
   {
      if( !SCIPisZero(scip, gcgcol->sepamastercutcoeffs[j]) )
         norm += SQR(gcgcol->sepamastercutcoeffs[j]);;
   }

   for( i = 0; i < gcgcol->ninferredpricingvars; ++i )
   {
      assert(!SCIPisZero(scip, gcgcol->inferredpricingvals[i]));
      norm += SQR(gcgcol->inferredpricingvals[i]);
   }

   for( i = 0; i < gcgcol->nlinkvars; ++i )
   {
      if( !SCIPisZero(scip, gcgcol->vals[gcgcol->linkvars[i]]) )
         norm += gcgcol->vals[gcgcol->linkvars[i]];
   }

   /* consider convexity constraint */
   norm += 1.0;

   gcgcol->norm = norm;
   return SCIP_OKAY;
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

/** update the column's stored coefficients for separator mastercut */
SCIP_RETCODE GCGcolAppendSepaMastercutCoeffs(
   GCG_COL*             gcgcol,                 /**< gcg column structure */
   SCIP_Real*           sepamastercoeffs,       /**< pointer to array of new mastercut coefficients */
   int                  nsepamastercoeffs       /**< number of new mastercut coefficients */
   )
{
   int i;

   SCIPdebugMessage("append %i coefficients for master separator cuts\n", nsepamastercoeffs);

   /* allocate memory for new coefficients (if necessary) */
   if( gcgcol->sepamastercutscoeffssize < gcgcol->nsepamastercutcoeffs + nsepamastercoeffs )
   {
      int newsize = SCIPcalcMemGrowSize(gcgcol->pricingprob, gcgcol->nsepamastercutcoeffs + nsepamastercoeffs);
      SCIP_CALL( SCIPreallocBlockMemoryArray(gcgcol->pricingprob, &(gcgcol->sepamastercutcoeffs), gcgcol->sepamastercutscoeffssize, newsize) );
      gcgcol->sepamastercutscoeffssize = newsize;
   }

   /* transfer coefficients to array */
   SCIPdebugMessage("current number of coefficients %i\n", gcgcol->nsepamastercutcoeffs);
   for( i = 0; i < nsepamastercoeffs; i++ )
   {
      SCIPdebugMessage("add %f at %i\n", sepamastercoeffs[i], gcgcol->nsepamastercutcoeffs);
      gcgcol->sepamastercutcoeffs[gcgcol->nsepamastercutcoeffs] = sepamastercoeffs[i];
      ++(gcgcol->nsepamastercutcoeffs);
   }

   return SCIP_OKAY;
}


/** get the column's number of stored coefficients for separator mastercuts */
int GCGcolGetNSepaMastercutCoeffs(
   GCG_COL*             gcgcol      /**< gcg column structure */
   )
{
   return gcgcol->nsepamastercutcoeffs;
}

/** get the size of the column's array for separator mastercut coefficients */
int GCGcolGetSepaMastercutCoeffsSize(
   GCG_COL*            gcgcol    /**< gcg column structure */
   )
{
   return gcgcol->sepamastercutscoeffssize;
}

/** get the column's coefficients for the separator mastercuts */
SCIP_Real* GCGcolGetSepaMastercutCoeffs(
   GCG_COL*             gcgcol      /**< gcg column structure */
   )
{
   return gcgcol->sepamastercutcoeffs;
}

/** return solution value of variable in gcg column */
SCIP_Real GCGcolGetSolVal(
   GCG_COL*             gcgcol,             /**< gcg column */
   SCIP_VAR*            var                 /**< variable */
   )
{
   SCIP_VAR** vars = NULL;
   SCIP_Real* vals;
   int nvars;
   int pos;

   if( GCGvarIsPricing(var) )
   {
      vars = gcgcol->vars;
      vals = gcgcol->vals;
      nvars = gcgcol->nvars;
   }
   else if( GCGvarIsInferredPricing(var) )
   {
      vars = gcgcol->inferredpricingvars;
      vals = gcgcol->inferredpricingvals;
      nvars = gcgcol->ninferredpricingvars;
   }

   if( vars != NULL && SCIPsortedvecFindPtr((void**) vars, SCIPvarComp, (void*) var, nvars, &pos) )
      return vals[pos];

   return 0.0;
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
SCIP_RETCODE GCGcolComputeDualObjPara(
   GCG*                 gcg,                /**< GCG data structure */
   GCG_COL*             gcgcol,             /**< gcg column */
   SCIP_Real*           para                /**< pointer to store the parallelism of column to dual objective */
   )
{
   int i;
   int j;

   GCG_EXTENDEDMASTERCONSDATA** activecuts;
   SCIP_CONS** masterconss;
   SCIP_ROW** originalsepamastercuts;
   SCIP_Real dualobjnorm;
   SCIP* masterprob;

   assert(gcg != NULL);
   
   masterprob = GCGgetMasterprob(gcg);

   assert(gcgcol != NULL);

   masterconss = GCGgetMasterConss(gcg);
   originalsepamastercuts = GCGsepaGetOriginalSepaMastercuts(gcg);

   *para = 0.0;

   dualobjnorm = 0.0;

   /* compute scalar of master values of gcg columns */
   for( i = 0; i < gcgcol->nmastercoefs; ++i )
   {
      SCIP_Real lhs;
      SCIP_Real rhs;

      lhs = SCIPgetLhsLinear(masterprob, masterconss[i]);
      rhs = SCIPgetRhsLinear(masterprob, masterconss[i]);

      if( !SCIPisInfinity(masterprob, -lhs) )
      {
         dualobjnorm += SQR(lhs);

         if( SCIPisPositive(masterprob, gcgcol->mastercoefs[i]) )
            *para += gcgcol->mastercoefs[i] * lhs;
      }
      else if( !SCIPisInfinity(masterprob, rhs) )
      {
         dualobjnorm += SQR(rhs);

         if(SCIPisNegative(masterprob, gcgcol->mastercoefs[i] ) )
            *para += gcgcol->mastercoefs[i] * rhs;
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

      if( !SCIPisInfinity(masterprob, -lhs))
      {
         dualobjnorm += SQR(lhs);

         if( SCIPisPositive(masterprob, gcgcol->originalsepamastercuts[i]) )
            *para += gcgcol->originalsepamastercuts[i] * lhs;
      }
      else if( !SCIPisInfinity(masterprob, rhs) )
      {
         dualobjnorm += SQR(rhs);

         if(SCIPisNegative(masterprob, gcgcol->originalsepamastercuts[i] ) )
            *para += gcgcol->originalsepamastercuts[i] * rhs;
      }
   }

   for( i = 0; i < gcgcol->ninferredpricingvars; ++i )
   {
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real coef;
      GCG_EXTENDEDMASTERCONSDATA* branchextendedmasterconsdata;

      coef = gcgcol->inferredpricingvals[i];
      assert(!SCIPisZero(masterprob, coef));

      branchextendedmasterconsdata = GCGinferredPricingVarGetExtendedmasterconsdata(gcgcol->inferredpricingvars[i]);
      lhs = GCGextendedmasterconsGetLhs(gcg, branchextendedmasterconsdata);
      rhs = GCGextendedmasterconsGetRhs(gcg, branchextendedmasterconsdata);

      if( !SCIPisInfinity(masterprob, -lhs))
      {
         dualobjnorm += SQR(lhs);

         if( SCIPisPositive(masterprob, coef) )
            *para += coef * lhs;
      }
      else if( !SCIPisInfinity(masterprob, rhs) )
      {
         dualobjnorm += SQR(rhs);

         if( SCIPisPositive(masterprob, coef) )
            *para += coef * rhs;
      }
   }

   /* consider cuts generated by separators */
   activecuts = GCGgetActiveCuts(gcg);
   for( j = 0; j < gcgcol->nsepamastercutcoeffs; j++ )
   {
      SCIP_ROW* mastercutrow;
      SCIP_Real lhs;
      SCIP_Real rhs;

      if( !GCGextendedmasterconsIsActive(activecuts[j]) )
         continue;

      mastercutrow = GCGextendedmasterconsGetRow(activecuts[j]);
      lhs = SCIProwGetLhs(mastercutrow);
      rhs = SCIProwGetRhs(mastercutrow);

      if( !SCIPisInfinity(masterprob, -lhs) )
      {
         dualobjnorm += SQR(lhs);

         if( SCIPisPositive(masterprob, gcgcol->sepamastercutcoeffs[j]) )
            *para += gcgcol->sepamastercutcoeffs[j] * lhs;
      }
      else if( !SCIPisInfinity(masterprob, rhs) )
      {
         dualobjnorm += SQR(rhs);

         if( SCIPisNegative(masterprob, gcgcol->sepamastercutcoeffs[j]) )
            *para += gcgcol->sepamastercutcoeffs[j] * rhs;
      }
   }

   for( i = 0; i < GCGgetNPricingprobs(gcg); ++i )
      dualobjnorm += SQR(GCGgetNIdenticalBlocks(gcg, i));

   *para += SQR(GCGgetNIdenticalBlocks(gcg, gcgcol->probnr));

   assert(!SCIPisInfinity(masterprob, ABS(*para)));

   dualobjnorm = sqrt(dualobjnorm);
   assert(!SCIPisInfinity(masterprob, dualobjnorm));
   assert(SCIPisPositive(masterprob, dualobjnorm));
   assert(SCIPisPositive(masterprob, gcgcol->norm));

   *para = (*para) / (dualobjnorm * gcgcol->norm);

   return SCIP_OKAY;
}

/** compute orthogonality of two gcg columns */
SCIP_RETCODE GCGcolComputeOrth(
   GCG*                 gcg,                /**< GCG data structure */
   GCG_COL*             gcgcol1,            /**< first gcg column */
   GCG_COL*             gcgcol2,            /**< second gcg column */
   SCIP_Real*           orth                /**< pointer to store the orthogonality of two gcg columns */
   )
{
   SCIP* scip;
   int i;
   int j;
   SCIP_Real para = 0.0;
   SCIP_Real norm1 = 0.0;
   SCIP_Real norm2 = 0.0;

   assert(gcg != NULL);
   assert(gcgcol1 != NULL);
   assert(gcgcol2 != NULL);

   scip = GCGgetMasterprob(gcg);

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

   /* consider cuts generated by separators */
   assert(gcgcol1->nsepamastercutcoeffs == gcgcol2->nsepamastercutcoeffs);
   for( j = 0; j < gcgcol1->nsepamastercutcoeffs; j++ )
   {
      if( SCIPisPositive(scip, gcgcol1->sepamastercutcoeffs[j] * gcgcol2->sepamastercutcoeffs[j]) )
         para += gcgcol1->sepamastercutcoeffs[j] * gcgcol2->sepamastercutcoeffs[j];

      if( SCIPisPositive(scip, gcgcol1->sepamastercutcoeffs[j]) )
         norm1 += SQR(gcgcol1->sepamastercutcoeffs[j]);

      if( SCIPisPositive(scip, gcgcol2->sepamastercutcoeffs[j]) )
         norm2 += SQR(gcgcol2->sepamastercutcoeffs[j]);
   }

   i = 0;
   j = 0;
   while( i < gcgcol1->ninferredpricingvars || j < gcgcol2->ninferredpricingvars )
   {
      if( i < gcgcol1->ninferredpricingvars && j < gcgcol2->ninferredpricingvars &&
          gcgcol1->inferredpricingvars[i] == gcgcol2->inferredpricingvars[j] )
      {
         SCIP_Real coef1 = gcgcol1->inferredpricingvals[i];
         SCIP_Real coef2 = gcgcol2->inferredpricingvals[j];

         assert(!SCIPisZero(scip, coef1));
         
         if( SCIPisPositive(scip, coef1 * coef2) )
            para += coef1 * coef2;

         if( SCIPisPositive(scip, coef1) )
            norm1 += SQR(coef1);
         if( SCIPisPositive(scip, coef2) )
            norm2 += SQR(coef2);
         ++i;
         ++j;
      }
      else if( i < gcgcol1->ninferredpricingvars && 
               SCIPvarGetIndex(gcgcol1->inferredpricingvars[i]) < SCIPvarGetIndex(gcgcol2->inferredpricingvars[j]) )
      {
         if( SCIPisPositive(scip, gcgcol1->inferredpricingvals[i]) )
            norm1 += SQR(gcgcol1->inferredpricingvals[i]);
         ++i;
      }
      else
      {
         assert(j < gcgcol2->ninferredpricingvars);
         assert(SCIPvarGetIndex(gcgcol1->inferredpricingvars[i]) > SCIPvarGetIndex(gcgcol2->inferredpricingvars[j]));
         if( SCIPisPositive(scip, gcgcol2->inferredpricingvals[j]) )
            norm2 += SQR(gcgcol2->inferredpricingvals[j]);
         ++j;
      }
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

   *orth = 1.0 - para;
   return SCIP_OKAY;
}

/** returns the inferred (coefficient) pricing variables solution values */
SCIP_Real* GCGcolGetInferredPricingVals(
   GCG_COL*              gcgcol             /**< gcgcol */
   )
{
   return gcgcol->inferredpricingvals;
}

/** returns the inferred (coefficient) pricing variables */
SCIP_VAR** GCGcolGetInferredPricingVars(
   GCG_COL*              gcgcol             /**< gcgcol */
   )
{
   return gcgcol->inferredpricingvars;
}

/** returns the number of inferred (coefficient) pricing variables */
int GCGcolGetNInferredPricingVars(
   GCG_COL*              gcgcol             /**< gcgcol */
   )
{
   return gcgcol->ninferredpricingvars;
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
