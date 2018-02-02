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

/**@file   gcgcol.c
 * @brief  methods for working with gcg column structure
 * @author Jonas Witt
 *
 * Various methods to work with gcg column structure
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "pub_gcgcol.h"

#include "gcg.h"
#include "scip/def.h"
#include "scip/scip.h"
#include "scip_misc.h"
#include "blockmemshell/memory.h"
#include "relax_gcg.h"

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

   SCIP_CALL( SCIPallocMemory(pricingprob, gcgcol) );


   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &((*gcgcol)->vars), nvars) );
   SCIP_CALL( SCIPallocMemoryArray(pricingprob, &((*gcgcol)->vals), nvars) );

   (*gcgcol)->pricingprob = pricingprob;
   (*gcgcol)->probnr = probnr;
   (*gcgcol)->isray = isray;
   (*gcgcol)->redcost = redcost;
   (*gcgcol)->age = 0;

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

      if( !SCIPisZero(pricingprob, origval) )
      {
         (*gcgcol)->vars[nnonz] = origvar;
         (*gcgcol)->vals[nnonz] = origval;
         ++nnonz;
      }
   }

   (*gcgcol)->nvars = nnonz;

   /* sort vars and vals array w.r.t. variable index */
   SCIPsortPtrReal((void**)(*gcgcol)->vars, (double*)(*gcgcol)->vals, SCIPvarComp, nnonz);

#ifndef NDEBUG
   for( i = 1 ; i < (*gcgcol)->nvars; ++i )
   {
      assert( SCIPvarCompare((*gcgcol)->vars[i-1], (*gcgcol)->vars[i]) != 0 );
   }
#endif
   return SCIP_OKAY;
}

/** free a gcg column */
void GCGfreeGcgCol(
   GCG_COL**            gcgcol              /**< pointer to store gcg column */
   )
{
   assert(gcgcol != NULL);
   assert(*gcgcol != NULL);

   /* todo: release vars? */
   assert((*gcgcol)->vars != NULL);
   SCIPfreeMemoryArray((*gcgcol)->pricingprob, &(*gcgcol)->vars);
   assert((*gcgcol)->vals != NULL);
   SCIPfreeMemoryArray((*gcgcol)->pricingprob, &(*gcgcol)->vals);
   SCIPfreeMemory((*gcgcol)->pricingprob, gcgcol);
}

/** create a gcg column from a solution to a pricing problem */
SCIP_RETCODE GCGcreateGcgColFromSol(
   SCIP*                pricingprob,        /**< SCIP data structure (pricing problem) */
   GCG_COL**            gcgcol,             /**< pointer to store gcg column */
   int                  prob,               /**< number of corresponding pricing problem */
   SCIP_SOL*            sol,                /**< solution of pricing problem with index prob */
   SCIP_Bool            isray,              /**< is column a ray? */
   SCIP_Real            redcost             /**< last known reduced cost */
)
{
   SCIP_VAR** solvars;
   SCIP_VAR** colvars;

   SCIP_Real* colvals;

   int nsolvars;
   int ncolvars;

   int i;

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
      solval = SCIPgetSolVal(pricingprob, sol, solvar);

      if( SCIPisZero(pricingprob, solval) )
      {
         continue;
      }

      colvars[ncolvars] = solvar;
      colvals[ncolvars] = solval;
      ++ncolvars;
   }

   SCIP_CALL( GCGcreateGcgCol(pricingprob, gcgcol, prob, colvars, colvals, ncolvars, isray, redcost) );

   SCIPfreeBufferArray(pricingprob, &colvars);
   SCIPfreeBufferArray(pricingprob, &colvals);

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
      SCIP_VAR* var1;
      SCIP_VAR* var2;

      SCIP_Real val1;
      SCIP_Real val2;

      var1 = vars1[i];
      var2 = vars2[i];

      val1 = vals1[i];
      val2 = vals2[i];

      if( SCIPvarCompare(var1, var2) != 0 || !SCIPisEQ(pricingprob, val1, val2) )
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
SCIP_RETCODE GCGcolUpdateRedcost(
   GCG_COL*             gcgcol,             /**< gcg column structure */
   SCIP_Real            redcost,            /**< new reduced cost */
   SCIP_Bool            growold             /**< change age depending on reduced cost? */
   )
{
   gcgcol->redcost = redcost;

   if( !growold )
   {
      return SCIP_OKAY;
   }

   if( !SCIPisNegative(gcgcol->pricingprob, redcost) )
   {
      ++(gcgcol->age);
   }
   else
   {
      gcgcol->age = 0;
   }

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
