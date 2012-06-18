/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"
//#define SCIP_DEBUG
/**@file   branch_generic.c
 * @ingroup BRANCHINGRULES
 * @brief  branching rule based on vanderbeck's generic branching scheme
 * @author Marcel Schmickerath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_generic.h"
#include "relax_gcg.h"
#include "cons_origbranch.h"
#include "pricer_gcg.h"
#include "scip/cons_varbound.h"
#include "type_branchgcg.h"
#include "struct_branchgcg.h"
#include "struct_vardata.h"

#include <stdio.h>
#include <stdlib.h>


#define BRANCHRULE_NAME          "generic"
#define BRANCHRULE_DESC          "generic branching rule by Vanderbeck"
#define BRANCHRULE_PRIORITY      10  //?
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

typedef int ComponentBoundSequence[3];  // [[comp], [sense], [bound]]


/** branching data for branching decisions */
struct GCG_BranchData
{
   SCIP_VAR*          var1;                  /**< first original variable on which the branching is done */
   SCIP_VAR*          var2;                  /**< second original variable on which the branching is done */
   SCIP_Bool          same;                  /**< should each master var contain either both or none of the vars? */
   int                blocknr;               /**< number of the block in which branching was performed */
   SCIP_CONS*         pricecons;             /**< constraint enforcing the branching restriction in the pricing problem */
};

struct GCG_BranchData
{
   SCIP_VAR*          mastervar;             /**< master variable */
   int                blocknr;               /**< number of the block in which branching was performed */
   SCIP_Real*         generator;             /**< corresponding generator to the mastervar */
   int                generatorsize;
   SCIP_CONS*         pricecons;             /**< constraint enforcing the branching restriction in the pricing problem */
};


/*
 * Callback methods for enforcing branching constraints
 */

static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterRyanfoster)
{
   SCIP* origscip;
   SCIP* pricingscip;
   SCIP_VARDATA* vardata1;
   SCIP_VARDATA* vardata2;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->var1 != NULL);
   assert(branchdata->var2 != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   pricingscip = GCGrelaxGetPricingprob(origscip, branchdata->blocknr);
   assert(pricingscip != NULL);

   SCIPdebugMessage("branchActiveMasterRyanfoster: %s(%s, %s)\n", ( branchdata->same ? "same" : "differ" ),
      SCIPvarGetName(branchdata->var1), SCIPvarGetName(branchdata->var2));

   /* get vardatas */
   vardata1 = SCIPvarGetData(branchdata->var1);
   assert(vardata1 != NULL);
   assert(vardata1->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata1->blocknr == branchdata->blocknr);
   assert(vardata1->data.origvardata.pricingvar != NULL);

   vardata2 = SCIPvarGetData(branchdata->var2);
   assert(vardata2 != NULL);
   assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata2->blocknr == branchdata->blocknr);
   assert(vardata2->data.origvardata.pricingvar != NULL);

   /* create corresponding constraint in the pricing problem, if not yet created */
   if( branchdata->pricecons == NULL )
   {
      if( branchdata->same )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "same(%s, %s)", SCIPvarGetName(branchdata->var1), 
            SCIPvarGetName(branchdata->var2));

         SCIP_CALL( SCIPcreateConsVarbound(pricingscip,
               &(branchdata->pricecons), name, vardata1->data.origvardata.pricingvar, 
               vardata2->data.origvardata.pricingvar, -1.0, 0.0, 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      }
      else
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "differ(%s, %s)", SCIPvarGetName(branchdata->var1), 
            SCIPvarGetName(branchdata->var2));

         SCIP_CALL( SCIPcreateConsVarbound(pricingscip,
               &(branchdata->pricecons), name, vardata1->data.origvardata.pricingvar, 
               vardata2->data.origvardata.pricingvar, 1.0, -SCIPinfinity(scip), 1.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      }
   }
   /* add constraint to the pricing problem that enforces the branching decision */
   SCIP_CALL( SCIPaddCons(pricingscip, branchdata->pricecons) );
   
   return SCIP_OKAY;
}

static
GCG_DECL_BRANCHDEACTIVEMASTER(branchDeactiveMasterRyanfoster)
{
   SCIP* origscip;
   SCIP* pricingscip;

   assert(scip != NULL);
    assert(branchdata != NULL);
   assert(branchdata->var1 != NULL);
   assert(branchdata->var2 != NULL);
   assert(branchdata->pricecons != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   pricingscip = GCGrelaxGetPricingprob(origscip, branchdata->blocknr);
   assert(pricingscip != NULL);

   SCIPdebugMessage("branchDeactiveMasterRyanfoster: %s(%s, %s)\n", ( branchdata->same ? "same" : "differ" ),
      SCIPvarGetName(branchdata->var1), SCIPvarGetName(branchdata->var2));

   /* remove constraint from the pricing problem that enforces the branching decision */
   assert(branchdata->pricecons != NULL);
   SCIP_CALL( SCIPdelCons(pricingscip, branchdata->pricecons) );
   
   return SCIP_OKAY;
}

static
GCG_DECL_BRANCHPROPMASTER(branchPropMasterRyanfoster)
{
   SCIP* origscip;
   SCIP_VARDATA* vardata;
   SCIP_VAR** vars;
   SCIP_Real val1;
   SCIP_Real val2;
   int nvars;
   int propcount;
   int i;
   int j;

   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->var1 != NULL);
   assert(branchdata->var2 != NULL);
   assert(branchdata->pricecons != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   SCIPdebugMessage("branchPropMasterRyanfoster: %s(%s, %s)\n", ( branchdata->same ? "same" : "differ" ),
      SCIPvarGetName(branchdata->var1), SCIPvarGetName(branchdata->var2));

   *result = SCIP_DIDNOTFIND;

   propcount = 0;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
          
   /* iterate over all master variables */
   for( i = 0; i < nvars; i++)
   {
      /* only look at variables not fixed to 0 */
      if( !SCIPisFeasZero(scip, SCIPvarGetUbLocal(vars[i])) )
      {
         vardata = SCIPvarGetData(vars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_MASTER);
         assert(vardata->blocknr >= -1 && vardata->blocknr < GCGrelaxGetNPricingprobs(origscip));
         assert(vardata->data.mastervardata.norigvars >= 0);
         assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
         assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);

         /* if variable belongs to a different block than the branching restriction, we do not have to look at it */
         if( branchdata->blocknr != vardata->blocknr )
            continue;
         
         /* save the values of the original variables for the current master variable */
         val1 = 0.0;
         val2 = 0.0;
         for( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
         {
            if( vardata->data.mastervardata.origvars[j] == branchdata->var1 )
            {
               assert(SCIPisEQ(scip, vardata->data.mastervardata.origvals[j], 1.0));
               val1 = vardata->data.mastervardata.origvals[j];
               continue;
            }
            if( vardata->data.mastervardata.origvars[j] == branchdata->var2 )
            {
               assert(SCIPisEQ(scip, vardata->data.mastervardata.origvals[j], 1.0));
               val2 = vardata->data.mastervardata.origvals[j];
            }
         }

         /* if branching enforces that both original vars are either both contained or none of them is contained
          * and the current master variable has different values for both of them, fix the variable to 0 */
         if( branchdata->same && !SCIPisEQ(scip, val1, val2) )
         {
            SCIPchgVarUb(scip, vars[i], 0.0);
            propcount++;
         }
         /* if branching enforces that both original vars must be in different mastervars, fix all 
          * master variables to 0 that contain both */
         if( !branchdata->same && SCIPisEQ(scip, val1, 1.0) && SCIPisEQ(scip, val1, 1.0) )
         {
            SCIPchgVarUb(scip, vars[i], 0.0);
            propcount++;
         }
      }
   }
      
   SCIPdebugMessage("Finished propagation of branching decision constraint: %s(%s, %s), %d vars fixed.\n", 
      ( branchdata->same ? "same" : "differ" ), SCIPvarGetName(branchdata->var1), SCIPvarGetName(branchdata->var2), propcount);

   if( propcount > 0 )
   {
      *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}

static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteRyanfoster)
{
   assert(scip != NULL);
   assert(branchdata != NULL);
   
   SCIPdebugMessage("branchDataDeleteRyanfoster: %s(%s, %s)\n", ( (*branchdata)->same ? "same" : "differ" ),
      SCIPvarGetName((*branchdata)->var1), SCIPvarGetName((*branchdata)->var2));

   /* release constraint that enforces the branching decision */
   if( (*branchdata)->pricecons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(GCGrelaxGetPricingprob(scip, (*branchdata)->blocknr), 
            &(*branchdata)->pricecons) );
   }

   SCIPfreeMemory(scip, branchdata);
   *branchdata = NULL;

   return SCIP_OKAY;
}

/*
 * Callback methods
 */




/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpRyanfoster)
{  
   SCIPdebugMessage("Execlp method of ryanfoster branching\n");
//   printf("Execlp method of ryanfoster branching\n");

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextRyanfoster)
{
   SCIP* masterscip;
   SCIP_VAR** mastervars;
   int nmastervars;
   int nbinmastervars;
   int nintmastervars;

   SCIP_NODE* childsame;
   SCIP_NODE* childdiffer;
   SCIP_CONS* origbranchsame;
   SCIP_CONS* origbranchdiffer;
   GCG_BRANCHDATA* branchsamedata;
   GCG_BRANCHDATA* branchdifferdata;
   char samename[SCIP_MAXSTRLEN];
   char differname[SCIP_MAXSTRLEN];

   SCIP_CONS* origcons;

   SCIP_Bool feasible;
   SCIP_Bool contained;

   SCIP_VAR** branchcands;
   SCIP_Real* branchcandsfrac;
   SCIP_Real* branchcandssol;
   int nbranchcands;

   SCIP_VARDATA* vardata1;
   SCIP_VARDATA* vardata2;

   int v1;
   int v2;
   int o1;
   int o2;
   int i;
   int j;

   SCIP_VAR* mvar1;
   SCIP_VAR* mvar2;
   SCIP_VAR* ovar1;
   SCIP_VAR* ovar2;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execrel method of ryanfoster branching\n");

   *result = SCIP_DIDNOTRUN;

   /* check whether the current original solution is integral */
   SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), TRUE, TRUE, TRUE, TRUE, &feasible) );
   if( feasible )
   {
      SCIPdebugMessage("node cut off, since origsol was feasible, solval = %f\n", SCIPgetSolOrigObj(scip, GCGrelaxGetCurrentOrigSol(scip)));
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* the current original solution is not integral, now we have to branch;
    * first, get the master problem and all variables of the master problem */
   masterscip = GCGrelaxGetMasterprob(scip);
   SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, &nbinmastervars, &nintmastervars, NULL, NULL) );
   SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, &branchcandssol, &branchcandsfrac, &nbranchcands, NULL) );

   /* now search for 2 (fractional) columns v1, v2 in the master and 2 original variables o1, o2
    * s.t. v1 contains both o1 and o2 and column 2 contains either o1 or o2
    */
   ovar1 = NULL;
   ovar2 = NULL;
   vardata1 = NULL;

   feasible = FALSE;
   for( v1 = 0; v1 < nbranchcands && !feasible; v1++ )
   {
      mvar1 = branchcands[v1];
      vardata1 = SCIPvarGetData(mvar1);
      assert(vardata1 != NULL);
      assert(vardata1->vartype == GCG_VARTYPE_MASTER);
      for( o1 = 0; o1 < vardata1->data.mastervardata.norigvars && !feasible; o1++ )
      {
         ovar1 = vardata1->data.mastervardata.origvars[o1];
         /* v1 contains o1, look for v2 */
         for( v2 = v1+1; v2 < nbranchcands && !feasible; v2++ )
         {
            mvar2 = branchcands[v2];
            vardata2 = SCIPvarGetData(mvar2);
            assert(vardata2 != NULL);
            assert(vardata2->vartype == GCG_VARTYPE_MASTER);
            contained = FALSE;
            for( j = 0; j < vardata2->data.mastervardata.norigvars; j++ )
            {
               if( vardata2->data.mastervardata.origvars[j] == ovar1 )
               {
                  contained = TRUE;
                  break;
               }
            }

            if(!contained)
               continue;

            /* v2 also contains o1, now look for o2 */
            for( o2 = 0; o2 < vardata1->data.mastervardata.norigvars && !feasible; o2++ )
            {
               ovar2 = vardata1->data.mastervardata.origvars[o2];
               if( ovar2 == ovar1 ) 
                  continue;

               contained = FALSE;
               for( j = 0; j < vardata2->data.mastervardata.norigvars; j++ )
               {
                  if( vardata2->data.mastervardata.origvars[j] == ovar2 )
                  {
                     contained = TRUE;
                     break;
                  }
               }

               if( contained )
                  continue;

               feasible = TRUE;
            }


            if( !feasible )
            {
               for( o2 = 0; o2 < vardata2->data.mastervardata.norigvars && !feasible; o2++ )
               {
                  ovar2 = vardata2->data.mastervardata.origvars[o2];
                  if( ovar2 == ovar1 ) 
                     continue;

                  contained = FALSE;
                  for( j = 0; j < vardata1->data.mastervardata.norigvars; j++ )
                  {
                     if( vardata1->data.mastervardata.origvars[j] == ovar2 )
                     {
                        contained = TRUE;
                        break;
                     }
                  }
                  
                  if( contained )
                     continue;
                  
                  feasible = TRUE;
               }
            }
         }
      }
   }

   if( !feasible )
   {
      SCIPdebugMessage("Ryanfoster branching rule could not find variables to branch on!\n");
      return SCIP_OKAY;
   }
   else
   {
      SCIPdebugMessage("Ryanfoster branching rule: branch on original variables %s and %s!\n",
         SCIPvarGetName(ovar1),
         SCIPvarGetName(ovar2));
   }

   assert(ovar1 != NULL);
   assert(ovar2 != NULL);
   assert(vardata1 != NULL);

   /* create the b&b-tree child-nodes of the current node */
   SCIP_CALL( SCIPcreateChild(scip, &childsame, 0.0, SCIPgetLocalTransEstimate(scip)) );
   SCIP_CALL( SCIPcreateChild(scip, &childdiffer, 0.0, SCIPgetLocalTransEstimate(scip)) );

   SCIP_CALL( SCIPallocMemory(scip, &branchsamedata) );
   SCIP_CALL( SCIPallocMemory(scip, &branchdifferdata) );

   branchsamedata->var1 = ovar1;
   branchsamedata->var2 = ovar2;
   branchsamedata->same = TRUE;
   branchsamedata->blocknr = vardata1->blocknr;
   branchsamedata->pricecons = NULL;

   branchdifferdata->var1 = ovar1;
   branchdifferdata->var2 = ovar2;
   branchdifferdata->same = FALSE;
   branchdifferdata->blocknr = vardata1->blocknr;
   branchdifferdata->pricecons = NULL;

   (void) SCIPsnprintf(samename, SCIP_MAXSTRLEN, "same(%s, %s)", SCIPvarGetName(branchsamedata->var1), 
      SCIPvarGetName(branchsamedata->var2));
   (void) SCIPsnprintf(differname, SCIP_MAXSTRLEN, "differ(%s, %s)", SCIPvarGetName(branchsamedata->var1),
      SCIPvarGetName(branchsamedata->var2));


   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchsame, samename, childsame, 
         GCGconsOrigbranchGetActiveCons(scip), branchrule, branchsamedata) );
   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchdiffer, differname, childdiffer, 
         GCGconsOrigbranchGetActiveCons(scip), branchrule, branchdifferdata) );

   /* add constraints to nodes */
   SCIP_CALL( SCIPaddConsNode(scip, childsame, origbranchsame, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childdiffer, origbranchdiffer, NULL) );

   /* add branching decision as linear constraints to original problem */
   vardata1 = SCIPvarGetData(branchdifferdata->var1);
   vardata2 = SCIPvarGetData(branchdifferdata->var2);
   assert(vardata1 != NULL);
   assert(vardata2 != NULL);
   assert(vardata1->vartype == GCG_VARTYPE_ORIGINAL);
   assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);

   vardata1 = SCIPvarGetData(vardata1->data.origvardata.pricingvar);
   vardata2 = SCIPvarGetData(vardata2->data.origvardata.pricingvar);
   assert(vardata1 != NULL);
   assert(vardata2 != NULL);
   assert(vardata1->vartype == GCG_VARTYPE_PRICING);
   assert(vardata2->vartype == GCG_VARTYPE_PRICING);
   assert(vardata1->blocknr == vardata2->blocknr);
   assert(vardata1->data.pricingvardata.norigvars == vardata2->data.pricingvardata.norigvars);

   for( i = 0; i < vardata1->data.pricingvardata.norigvars; i++ )
   {
      assert(SCIPvarGetData(vardata1->data.pricingvardata.origvars[i])->blocknr 
         == SCIPvarGetData(vardata2->data.pricingvardata.origvars[i])->blocknr);

      /* create constraint for same-child */
      SCIP_CALL( SCIPcreateConsVarbound(scip, &origcons, samename, 
            vardata1->data.pricingvardata.origvars[i], vardata2->data.pricingvardata.origvars[i], 
            -1.0, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      
      /* add cons locally to the problem */
      SCIP_CALL( SCIPaddConsNode(scip, childsame, origcons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &origcons) );
      assert(origcons == NULL);

      /* create constraint for differ-child */
      SCIP_CALL( SCIPcreateConsVarbound(scip, &origcons, differname, 
            vardata1->data.pricingvardata.origvars[i], vardata2->data.pricingvardata.origvars[i],
            1.0, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddConsNode(scip, childdiffer, origcons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &origcons) );
      assert(origcons == NULL);
   }

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &origbranchsame) );
   SCIP_CALL( SCIPreleaseCons(scip, &origbranchdiffer) );
      
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsRyanfoster)
{  
   SCIPdebugMessage("Execps method of ryanfoster branching\n");
   assert(0);

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitRyanfoster)
{  
   assert(branchrule != NULL);

   SCIP_CALL( GCGrelaxIncludeBranchrule(scip, branchrule, branchActiveMasterRyanfoster, 
         branchDeactiveMasterRyanfoster, branchPropMasterRyanfoster, NULL, branchDataDeleteRyanfoster) );
   
   return SCIP_OKAY;
}



/* define not used callback as NULL*/
#define branchCopyRyanfoster NULL
#define branchFreeRyanfoster NULL
#define branchExitRyanfoster NULL
#define branchInitsolRyanfoster NULL
#define branchExitsolRyanfoster NULL


/*
 * branching specific interface methods
 */

/** method for calculating the median over all fractional components values */
static
SCIP_Real GetMedian(SCIP_Real* array, int arraysize)
{  
	SCIP_Real Median;
	SCIP_Real swap;
	int l;
	int r;
	int i;
	int j;
	int MedianIndex;
	
	r = arraysize -1;
	l = 0;
	
	if( arraysize & 1)
		MedianIndex = arraysize/2;
	else
		MedianIndex = arraysize/2 -1;
	
	while(l < r-1)
	{
		Median = array[ MedianIndex ];
		i = l;
		j = r;
		do
		{
			while( array[i] < Median )
				++i;
			while( array[j] > Median )
				--j;
			if( i <=j )
			{
				swap = array[i];
			    array[i] = array[j];
			    array[j] = swap;
			    ++i;
			    --j;
			}
		} while( i <=j );
		if( j < MedianIndex )
			l = i;
		if( i > MedianIndex )
			r = j;
	}
	Median = array[ MedianIndex ];
	
   return Median;
}

// comparefunction for lexicographical sort
static
SCIP_DECL_SORTPTRCOMP(ptrcomp)
{
   struct GCG_BranchData* bd1;
   struct GCG_BranchData* bd2;
   int i;
   
   bd1 = (struct GCG_BranchData*) elem1;
   bd2 = (struct GCG_BranchData*) elem2;

   i = 0;
   
   for( i=0; i< bd1->generatorsize; ++i)
   {
	   if( bd1->generator[i] > bd2->generator[i] )
		   return -1;
	   if( bd1->generator[i] < bd2->generator[i] )
	   		   return 1;
   }
     
   return 0;
}

// lexicographical sort using qsort
static
SCIP_RETCODE LexicographicSort( GCG_BranchData* array, int arraysize)
{
   qsort( array, arraysize, sizeof(GCG_BranchData), ptrcomp );
   
   return SCIP_OKAY;
}


// compare function for ILO: returns true if bd1 < bd2 
static
SCIP_Bool ILOcomp( GCG_BranchData* bd1, GCG_BranchData* bd2, ComponentBoundSequence** C, int NBoundsequences, int* sequencesizes, ComponentBoundSequence* S, int Ssize, int* IndexSet, int indexsetsize, int p)
{
	int i;
	int j;
	int k;
	int l;
	int medianvalue;
	int newCsize;
	SCIP_Bool inall;
	SCIP_Bool inCj;
	SCIP_Bool inI;
	SCIP_Bool returnvalue;
	ComponentBoundSequence newcompbound;
	ComponentBoundSequence** copyC;
	int* copyI;
	
	i = -1;
	j = 0;
	k = 0;
	l = -1;
	inall = FALSE;
	inCj = FALSE;
	inI = FALSE;
	newCsize = 0;
	
	//lexicographic Order ?
	if( C == NULL )
		return (*ptrcomp( bd1, bd2) == -1);
  
	//find i which is in all S in C
	while(!inall)
	{
		++l;
		assert( l< indexsetsize);
		i = IndexSet[l];
	/*	inI = FALSE;
		for(j=0; j<indexsetsize; ++j)
		{
			if(IndexSet[j] == i)
			{
				inI =TRUE;
				break;
			}
		}
		if(!inI)
			continue;
		*/
		inall = TRUE;
		for(j=0; j< NBoundsequences; ++j)
		{
			inCj = FALSE;
			for(k=0; k<sequencesizes[j]; ++k)
			{
				if(C[j][k][0] == i)
				{
					inCj = TRUE;
					medianvalue = C[j][k][2];
					break;
				}
			}
			if(!inCj)
			{
				inall = FALSE;
				break;
			}
		}
	}
   	
	assert(i>=0);
	assert(i<indexsetsize);
	
	if(Ssize == 0)
		SCIP_CALL( SCIPallocBufferArray(scip, &S, 1) );
	else{
		//realloc S
		SCIP_CALL( SCIPallocBufferArray(scip, &copyS, Ssize) );
		
		for(j=0; j<Ssize; ++j)
			copyS[j]=S[j];
		
		SCIP_CALL( SCIPreallocBufferArray(scip, &S, Ssize+1) );
		
		for(j=0; j<Ssize; ++j)
			S[j]=copyS[j];
		
		SCIPfreeBufferArray(scip, &copyS);
	}
	++Ssize;
	
	//realloc I
	SCIP_CALL( SCIPallocBufferArray(scip, &copyI, indexsetsize) );

	for(j=0; j<indexsetsize; ++j)
		copyI[j]=IndexSet[j];

	SCIP_CALL( SCIPreallocBufferArray(scip, &IndexSet, indexsetsize-1) );

	k = 0;
	for(j=0; j<indexsetsize; ++j)
	{
		if(copyI[j] != i )
		{
			IndexSet[k]=copyI[j];
			++k;
		}
	}
	--indexsetsize;
	SCIPfreeBufferArray(scip, &copyI);

	
	//calculate subset of C
		   for(j=0; j< NBoundsequences; ++j)
		   {
			   for(k=0; k<sequencesizes[j]; ++k)
			   {
				   if(C[j][k][0] == i)
				   {
					   if(C[j][k[1] == 1)
						   ++Nupper;
					   else 
						   ++Nlower;
					   break;
				   }
			   }
		   }
		
   if( bd1->generator[i] == 1 && bd2->generator == 1 )
   {
	   newcompbound[0] = i;
	   newcompbound[1] = 1;
	   newcompbound[2] = medianvalue;
	   S[Ssize]=newbound;
	   k=0;
	   SCIP_CALL( SCIPreallocBufferArray(scip, &copyC, Nupper) );
	   SCIP_CALL( SCIPreallocBufferArray(scip, &newsequencesizes, Nupper) );
	   for(j=0; j< NBoundsequences; ++j)
	   {

		   for(k=0; k<sequencesizes[j]; ++k)
		   {
			   if(C[j][k][0] == i)
			   {
				   if(C[j][k][1] == 1)
				   {
					   copyC[k]=C[j];
					   newsequencesizes[k]=sequencesizes[j];
					   ++k;
					   break;
				   }
			   }
		   }
	   }

	   SCIP_CALL( SCIPreallocBufferArray(scip, &C, Nupper) );
	   
	   for(j=0;j<Nupper;++j)
		   C[j]=copyC[j];
	   
	   SCIPfreeBufferArray(scip, &copyC);

	   returnvalue = ILOcomp( bd1, bd2, C, Nupper, newsequencesizes, S, Ssize, IndexSet, indexsetsize, p+1);
   
	   SCIPfreeBufferArray(scip, &newsequencesizes);
	   
	   return returnvalue;
   }
   if( bd1->generator[i] == 0 && bd2->generator == 0 )
      {
	   newcompbound[0] = i;
	   newcompbound[1] = 0;
	   newcompbound[2] = medianvalue;
	   S[Ssize]=newbound;
	   k=0;
	   	   SCIP_CALL( SCIPreallocBufferArray(scip, &copyC, Nlower) );
	   	   SCIP_CALL( SCIPreallocBufferArray(scip, &newsequencesizes, Nlower) );
	   	   for(j=0; j< NBoundsequences; ++j)
	   	   {

	   		   for(k=0; k<sequencesizes[j]; ++k)
	   		   {
	   			   if(C[j][k][0] == i)
	   			   {
	   				   if(C[j][k][1] == 0)
	   				   {
	   					   copyC[k]=C[j];
	   					   newsequencesizes[k]=sequencesizes[j];
	   					   ++k;
	   					   break;
	   				   }
	   			   }
	   		   }
	   	   }

	   	   SCIP_CALL( SCIPreallocBufferArray(scip, &C, Nlower) );
	   	   
	   	   for(j=0;j<Nlower;++j)
	   		   C[j]=copyC[j];
	   	   
	   	   SCIPfreeBufferArray(scip, &copyC);

	   	   returnvalue = ILOcomp( bd1, bd2, C, Nlower, newsequencesizes, S, Ssize, IndexSet, indexsetsize, p+1);
	      
	   	   SCIPfreeBufferArray(scip, &newsequencesizes);
	   	   
	   	   return returnvalue; 
      }
   return (bd1->generator[i] > bd2->generator[i]);
}

// induced lexicographical sort based on QuickSort
static
SCIP_RETCODE ILOQSort( SCIP* scip, GCG_BranchData** array, int arraysize, ComponentBoundSequence** C, int sequencesize, int l, int r )
{
	int i;
	int j;
	int k;
	GCG_Branchdata* pivot;
	GCG_Branchdata* swap;
	int* IndexSet;
	indexsetsize;
	
	i = l;
	j = r;
	pivot = array[(l+r)/2];
	indexsetsize = pivot->generatorsize;
	k = 0;
	
	SCIP_CALL( SCIPallocBufferArray(scip, &IndexSet, indexsetsize) );
	for( k = 0; k < indexsetsize; ++k )
		IndexSet[k] = k; // ! n-1 here, instead of n
	
	do
	{
		while( ILOcomp( array[i], pivot, C, sequencesize, NULL, IndexSet, indexsetsize, 1))
			++i;
		while( ILOcomp( pivot, array[j], C, sequencesize, NULL, IndexSet, indexsetsize, 1))
					--j;
		if( i <= j )
		{
			swap = array[i];
			array[i] = array[j];
			array[j] = swap;
			++i;
			--j;
		}
	}while( i <= j );
	if( l < j )
		ILOQSort( scip, array, arraysize, C, sequencesize, l, j );
	if( i < r )
		ILOQSort( scip, array, arraysize, C, sequencesize, i, r );
	
	
	SCIPfreeBufferArray(scip, &IndexSet);
   
   return SCIP_OKAY;
}


// induced lexicographical sort
static
SCIP_RETCODE InducedLexicographicSort( SCIP* scip, GCG_BranchData* array, int arraysize, ComponentBoundSequence* C, int sequencesize )
{
	if( sequencesize == 0 )
		return LexicographicSort( array, arraysize );
	assert( S!= NULL );
	
   ILOQSort( scip, array, arraysize, C, sequencesize, 0, arraysize-1 );
   
   return SCIP_OKAY;
}

/** creates the most infeasible LP braching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleRyanfoster(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create branching rule data */
   branchruledata = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchCopyRyanfoster,
         branchFreeRyanfoster, branchInitRyanfoster, branchExitRyanfoster, branchInitsolRyanfoster, 
         branchExitsolRyanfoster, branchExeclpRyanfoster, branchExecextRyanfoster, branchExecpsRyanfoster,
         branchruledata) );

   return SCIP_OKAY;
}
