/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident ""

/**@file   heur_greedycolsel.c
 * @ingroup PRIMALHEURISTICS
 * @brief  greedy column selection primal heuristic
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* toggle debug mode */
//#define SCIP_DEBUG

#include <assert.h>

#include "heur_greedycolsel.h"
#include "pricer_gcg.h"
#include "pub_gcgvar.h"
#include "relax_gcg.h"


#define HEUR_NAME             "greedycolsel"
#define HEUR_DESC             "greedy column selection heuristic"
#define HEUR_DISPCHAR         'e'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             2
#define HEUR_FREQOFS          1
#define HEUR_MAXDEPTH         -1
/* TODO: should heuristic be called during the pricing loop or only after solving a node relaxation? */
#define HEUR_TIMING           SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_DURINGPRICINGLOOP
#define HEUR_USESSUBSCIP      FALSE

#define DEFAULT_MINCOLUMNS    200             /**< minimum number of columns to regard in the master problem */




/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                  mincolumns;           /**< minimum number of columns to regard in the master problem */

   int                  lastncols;            /**< number of columns in the last call of the heuristic       */
};




/*
 * Local methods
 */


/* How would the number of violated rows change if mastervar were increased?  */
static
int getViolationChange(
      SCIP*                   scip,
      SCIP_Real*              activities,
      SCIP_VAR*               mastervar
      )
{
   SCIP_COL* col;
   SCIP_ROW** colrows;
   SCIP_Real* colvals;
   int ncolrows;
   int violchange;

   int r;

   /* get the rows in which the master variable appears (only these must be regarded) */
   col = SCIPvarGetCol(mastervar);
   colrows = SCIPcolGetRows(col);
   colvals = SCIPcolGetVals(col);
   ncolrows = SCIPcolGetNLPNonz(col);
   assert(ncolrows == 0 || (colrows != NULL && colvals != NULL));

   violchange = 0;
   for( r = 0; r < ncolrows; r++ )
   {
      SCIP_ROW* row;
      int rowpos;

      row = colrows[r];
      rowpos = SCIProwGetLPPos(row);
      assert(-1 <= rowpos);

      if( rowpos >= 0 && !SCIProwIsLocal(row) )
      {
         SCIP_Real oldactivity;
         SCIP_Real newactivity;

         oldactivity = activities[rowpos];
         newactivity = oldactivity + colvals[r];

         if( SCIPisFeasLT(scip, oldactivity, SCIProwGetLhs(row)) || SCIPisFeasGT(scip, oldactivity, SCIProwGetRhs(row)) )
         {
            if( SCIPisFeasGE(scip, newactivity, SCIProwGetLhs(row)) && SCIPisFeasLE(scip, oldactivity, SCIProwGetRhs(row)) )
               violchange--;
         }
         else
         {
            if( SCIPisFeasLT(scip, newactivity, SCIProwGetLhs(row)) || SCIPisFeasGT(scip, newactivity, SCIProwGetRhs(row)) )
               violchange++;
         }
      }
   }

   return violchange;
}

/* get the index of the "best" master variable w.r.t. pseudo costs */
static
SCIP_RETCODE getBestMastervar(
      SCIP*                   scip,
      SCIP_Real*              activities,
      int*                    blocknr,
      int*                    index,
      int*                    violchange
      )
{
   SCIP* origprob;
   SCIP_VAR** mastervars;
   int nmastervars;

   SCIP_VAR* mastervar;
   int block;

   int i;
//   int j;
   int tmpviolchange;

   /* get original problem */
   origprob = GCGpricerGetOrigprob(scip);
   assert( origprob != NULL );

   /* get variable data of the master problem */
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(nmastervars >= 0);

   *index = -1;
   *violchange= SCIPinfinity(scip);

//   j = nmastervars - 1;
//   do
//   {
//      *index = j;
//      *violchange = getViolationChange(scip, activities, mastervars[j]);
//
//      vardata = SCIPvarGetData(mastervars[j]);
//      assert(vardata != NULL);
//      assert(vardata->vartype == GCG_VARTYPE_MASTER);
//
//      j--;
//   }
//   while( blocknr[vardata->blocknr] >= GCGrelaxGetNIdenticalBlocks(origprob, vardata->blocknr) );

   for( i = nmastervars - 1; i >= 0; i-- )
   {
      mastervar = mastervars[i];
      assert(GCGvarIsMaster(mastervar));
      block = GCGvarGetBlock(mastervar);

      /* TODO: handle copied original variables and linking variables */
      if( block < 0 )
         continue;

      /* ignore the master variable if the corresponding block is already full */
      if( blocknr[block] < GCGrelaxGetNIdenticalBlocks(origprob, block)
            && !GCGmasterVarIsRay(mastervar) ) /* TODO: handle rays */
      {
         tmpviolchange = getViolationChange(scip, activities, mastervar);
         if( tmpviolchange < *violchange )
         {
            *index = i;
            *violchange = tmpviolchange;
         }
      }
   }

   return SCIP_OKAY;
}

/* update activities */
static
SCIP_RETCODE updateActivities(
      SCIP*                   scip,
      SCIP_Real*              activities,
      SCIP_VAR*               mastervar
      )
{
   SCIP_COL* col;
   SCIP_ROW** colrows;
   SCIP_Real* colvals;
   int ncolrows;

   int r;

   assert(activities != NULL);

   col = SCIPvarGetCol(mastervar);
   colrows = SCIPcolGetRows(col);
   colvals = SCIPcolGetVals(col);
   ncolrows = SCIPcolGetNLPNonz(col);
   assert(ncolrows == 0 || (colrows != NULL && colvals != NULL));

   for( r = 0; r < ncolrows; ++r )
   {
      SCIP_ROW* row;
      int rowpos;

      row = colrows[r];
      rowpos = SCIProwGetLPPos(row);
      assert(-1 <= rowpos);

      if( rowpos >= 0 && !SCIProwIsLocal(row) )
      {
         SCIP_Real oldactivity;
         SCIP_Real newactivity;

         assert(SCIProwIsInLP(row));

         /* update row activity */
         oldactivity = activities[rowpos];
         newactivity = oldactivity + colvals[r];
         if( SCIPisInfinity(scip, newactivity) )
            newactivity = SCIPinfinity(scip);
         else if( SCIPisInfinity(scip, -newactivity) )
            newactivity = -SCIPinfinity(scip);
         activities[rowpos] = newactivity;
      }
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of primal heuristic
 */


/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#define heurCopyGreedycolsel NULL

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeGreedycolsel)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitGreedycolsel)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   heurdata->lastncols = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitGreedycolsel NULL


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolGreedycolsel NULL


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolGreedycolsel NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecGreedycolsel)
{  /*lint --e{715}*/
   SCIP* origprob;                           /* SCIP structure of original problem  */
   SCIP_HEURDATA* heurdata;                  /* heuristic's data                    */
   SCIP_ROW** lprows;                        /* LP rows of master problem           */
   SCIP_SOL* origsol;                        /* working original solution           */
   SCIP_VAR** mastervars;
   SCIP_VAR** origvars;
   SCIP_VAR** origpricingvars;
   SCIP_VAR* mastervar;
   SCIP_VAR* pricingvar;
   SCIP_Real* activities;                    /* for each master LP row, activity of current master solution          */
   SCIP_Real* origvals;
   int* blocknr;                             /* for each pricing problem, block we are currently working in          */
   SCIP_Bool allblocksfull;                  /* indicates if all blocks are full, i.e. all convexity constraints are satisfied */
   SCIP_Bool discretization;
   SCIP_Bool success;
   int block;
   int minnewcols;                           /* minimum number of new columns necessary for calling the heuristic    */
   int nlprows;
   int nmastervars;
   int norigvars;
   int norigpricingvars;
   int npricingprobs;
   int nviolrows;
   int violchange;

   int i;
   int index;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   /* get original problem */
   origprob = GCGpricerGetOrigprob(scip);
   assert( origprob != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   *result = SCIP_DIDNOTRUN;

   /* this heuristic works only for the discretization approach */
   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if( !discretization )
      return SCIP_OKAY;

   *result = SCIP_DELAYED;

   /* get variable data of the master problem */
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(nmastervars >= 0);

   /* calculate minimum number of new columns necessary for calling the heuristic;
    * this number is influenced by how successful the heuristic was in the past */
   minnewcols = heurdata->mincolumns * (int) 1.0 * ((1.0 + SCIPheurGetNCalls(heur)) / (1.0 + SCIPheurGetNBestSolsFound(heur)));

   /* if there are not enough new columns since last call, abort heuristic */
   if( nmastervars - heurdata->lastncols < minnewcols )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("Executing GCG greedy column selection heuristic (nmastervars = %d) ...\n", nmastervars);

   /* get number of pricing problems */
   npricingprobs = GCGrelaxGetNPricingprobs(origprob);
   assert( npricingprobs >= 0 );

   /* initialize the block numbers for the pricing problems */
   SCIP_CALL( SCIPallocBufferArray(scip, &blocknr, npricingprobs) );
   for( i = 0; i < npricingprobs; i++ )
   {
      blocknr[i] = 0;
   }
   allblocksfull = FALSE;

   /* get master LP rows data */
   SCIP_CALL( SCIPgetLPRowsData(scip, &lprows, &nlprows) );
   assert( lprows != NULL );
   assert( nlprows >= 0);

   /* get memory for working original solution and row activities */
   SCIP_CALL( SCIPcreateSol(origprob, &origsol, heur) );
   SCIP_CALL( SCIPallocBufferArray(scip, &activities, nlprows) );

   /* initialize activities with zero and get number of violated rows of zero master solution */
   nviolrows = 0;
   for( i = 0; i < nlprows; i++ )
   {
      SCIP_ROW* row;

      row = lprows[i];
      assert(SCIProwGetLPPos(row) == i);

      if( !SCIProwIsLocal(row) )
      {
         activities[i] = 0;
         if( SCIPisFeasLT(scip, 0, SCIProwGetLhs(row)) || SCIPisFeasGT(scip, 0, SCIProwGetRhs(row)) )
            nviolrows++;
      }
   }

   success = FALSE;

   /* try to increase master variables until all blocks are full */
   while( !allblocksfull && !success )
   {
      SCIP_CALL( getBestMastervar(scip, activities, blocknr, &index, &violchange) );

      /* if no master variable could be selected, abort */
      if( index == -1 )
         break;

      /* get master variable */
      mastervar = mastervars[index];
      assert(GCGvarIsMaster(mastervar));
      assert(GCGmasterVarIsRay(mastervar));

      /* get blocknr and original variables */
      block = GCGvarGetBlock(mastervar);
      origvars = GCGmasterVarGetOrigvars(mastervar);
      origvals = GCGmasterVarGetOrigvals(mastervar);
      norigvars = GCGmasterVarGetNOrigvars(mastervar);

      /* increase master value by one, i.e. increase solution values in current original solution accordingly */
      /* TODO: handle copied original variables and linking variables */
      if( block == -1 )
      {
         assert(norigvars == 1);
         assert(origvals[0] == 1.0);

         /* increase the corresponding value */
         SCIP_CALL( SCIPincSolVal(origprob, origsol, origvars[0], origvals[0]) );

         /* try to add original solution to solution pool */
         SCIP_CALL( SCIPtrySol(origprob, origsol, FALSE, TRUE, TRUE, TRUE, &success) );
      }
      else
      {
         /* loop over all original variables contained in the current master variable */
         for( i = 0; i < norigvars; i++ )
         {
            assert(!SCIPisZero(scip, origvals[i]));
            assert(GCGvarIsOriginal(origvars[i]));

            if(GCGvarGetBlock(origvars[i]) == -2)
               continue;

            pricingvar = GCGoriginalVarGetPricingVar(origvars[i]);
            assert(GCGvarIsPricing(pricingvar));

            norigpricingvars = GCGpricingVarGetNOrigvars(pricingvar);
            origpricingvars = GCGpricingVarGetOrigvars(pricingvar);

            /* increase the corresponding value */
            SCIP_CALL( SCIPincSolVal(origprob, origsol, origpricingvars[blocknr[block]], origvals[i]) );
         }

         blocknr[block]++;

         /* try to add original solution to solution pool */
         SCIP_CALL( SCIPtrySol(origprob, origsol, FALSE, TRUE, TRUE, TRUE, &success) );
      }

      /* update number of violated rows and activities array */
      nviolrows += violchange;
      SCIP_CALL( updateActivities(scip, activities, mastervars[index]) );

      /* check if all blocks are full */
      allblocksfull = TRUE;
      for( i = 0; i < npricingprobs; i++ )
      {
         allblocksfull &= blocknr[i] >= GCGrelaxGetNIdenticalBlocks(origprob, i);
      }
   }

//#ifdef SCIP_DEBUG
//   SCIPdebugMessage("  -> generated solution:\n");
//   SCIPprintSol(origprob, origsol, NULL, FALSE);
//#endif

   if( success )
   {
      *result = SCIP_FOUNDSOL;
      SCIPdebugMessage("  -> heuristic successful - feasible solution found.\n");
   }
   else
   {
      SCIPdebugMessage("  -> no feasible solution found.\n");
   }

   SCIPfreeSol(origprob, &origsol);
   SCIPfreeBufferArray(scip, &activities);
   SCIPfreeBufferArray(scip, &blocknr);

   heurdata->lastncols = nmastervars;

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the greedy column selection primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurGreedycolsel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create greedy column selection primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyGreedycolsel, heurFreeGreedycolsel, heurInitGreedycolsel, heurExitGreedycolsel,
         heurInitsolGreedycolsel, heurExitsolGreedycolsel, heurExecGreedycolsel,
         heurdata) );

   /* add greedy column selection primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/greedycolsel/mincolumns",
         "minimum number of columns to regard in the master problem",
         &heurdata->mincolumns, FALSE, DEFAULT_MINCOLUMNS, 1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
