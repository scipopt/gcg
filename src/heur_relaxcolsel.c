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

/**@file   heur_relaxcolsel.c
 * @ingroup PRIMALHEURISTICS
 * @brief  relaxation based column selection primal heuristic
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* toggle debug mode */
//#define SCIP_DEBUG

#include <assert.h>

#include "heur_relaxcolsel.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"
#include "struct_vardata.h"


#define HEUR_NAME             "relaxcolsel"
#define HEUR_DESC             "column selection heuristic that tries to round a master LP solution in promising directions"
#define HEUR_DISPCHAR         'x'
#define HEUR_PRIORITY         -100
#define HEUR_FREQ             2
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
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

/* initialize current working original solution as transformation of rounded down master LP solution
 * and add master variable candidates for rounding up to a list */
static
SCIP_RETCODE initializeOrigsol(
      SCIP*                   scip,
      SCIP_HEUR*              heur,
      SCIP_SOL*               origsol,
      SCIP_VAR***             mastercands,
      int*                    nmastercands,
      SCIP_Real*              activities,
      int*                    blocknr,
      SCIP_Bool*              success
      )
{
   SCIP* origprob;
   SCIP_VAR** mastervars;
   SCIP_VARDATA* vardata;
   SCIP_VARDATA* vardata2;
   SCIP_Real* mastervals;
   int nmastervars;
   int npricingprobs;

   int i;
   int j;

   /* get original problem */
   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   /* get variable data of the master problem */
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(nmastervars >= 0);

   /* get number of pricing problems */
   npricingprobs = GCGrelaxGetNPricingprobs(origprob);
   assert( npricingprobs >= 0 );

   SCIP_CALL( SCIPallocBufferArray(scip, &mastervals, nmastervars) );
   SCIP_CALL( SCIPgetSolVals(scip, NULL, nmastervars, mastervars, mastervals) );

   /* loop over all given master variables */
   for( i = 0; i < nmastervars; i++ )
   {
      vardata = SCIPvarGetData(mastervars[i]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_MASTER);
      assert(vardata->data.mastervardata.norigvars >= 0);
      assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
      assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);

      /* first of all, handle variables representing rays */
      if( vardata->data.mastervardata.isray )
      {
         assert(vardata->blocknr >= 0);
         /* we also want to take into account variables representing rays, that have a small value (between normal and feas eps),
          * so we do no feas comparison here */
         if( SCIPisPositive(scip, mastervals[i]) )
         {
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
            {
               assert(!SCIPisZero(scip, vardata->data.mastervardata.origvals[j]));

               assert(SCIPvarGetData(vardata->data.mastervardata.origvars[j]) != NULL);
               assert(SCIPvarGetData(vardata->data.mastervardata.origvars[j])->blocknr >= -2);
               assert(SCIPvarGetData(vardata->data.mastervardata.origvars[j])->blocknr < npricingprobs);

               /* the original variable is a linking variable */
               if( SCIPvarGetData(vardata->data.mastervardata.origvars[j])->blocknr == -2 )
                  continue;

               /* increase the corresponding value */
               SCIP_CALL( SCIPincSolVal(origprob, origsol, vardata->data.mastervardata.origvars[j], vardata->data.mastervardata.origvals[j] * SCIPfeasFloor(scip, mastervals[i])) );
            }
         }
         mastervals[i] = 0.0;
         continue;
      }

      /* handle the variables with value >= 1 to get integral values in original solution */
      /* TODO: handle copied original variables and linking variables */
      while( SCIPisFeasGE(scip, mastervals[i], 1) )
      {
         if( vardata->blocknr == -1 )
         {
            assert(vardata->data.mastervardata.norigvars == 1);
            assert(vardata->data.mastervardata.origvals[0] == 1.0);

            /* increase the corresponding value */
            SCIP_CALL( SCIPincSolVal(origprob, origsol, vardata->data.mastervardata.origvars[0], vardata->data.mastervardata.origvals[0]) );
            mastervals[i] = mastervals[i] - 1.0;
         }
         else
         {
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
            {
               assert(!SCIPisZero(scip, vardata->data.mastervardata.origvals[j]));

               /* get the right original variable */
               vardata2 = SCIPvarGetData(vardata->data.mastervardata.origvars[j]);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);

               if(vardata2->blocknr == -2)
                  continue;

               assert(vardata2->data.origvardata.pricingvar != NULL);
               vardata2 = SCIPvarGetData(vardata2->data.origvardata.pricingvar);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_PRICING);

               /* increase the corresponding value */
               SCIP_CALL( SCIPincSolVal(origprob, origsol, vardata2->data.pricingvardata.origvars[blocknr[vardata->blocknr]], vardata->data.mastervardata.origvals[j]) );
            }
            mastervals[i] = mastervals[i] - 1.0;
            blocknr[vardata->blocknr]++;
         }

         SCIP_CALL( updateActivities(scip, activities, mastervars[i]) );
      }

      /* if there is a fractional value >= 0.5 remaining for the master variable, add it as a candidate for rounding up */
      if( SCIPisFeasGE(scip, mastervals[i], 0.5)
            && !vardata->data.mastervardata.isray /* TODO: handle rays */
            && vardata->blocknr >= 0) /* TODO: handle copied original variables and linking variables */
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, mastercands, *nmastercands + 1) );
         (*mastercands)[*nmastercands] = mastervars[i];
         (*nmastercands)++;
      }
   }

   SCIP_CALL( SCIPtrySol(origprob, origsol, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &mastervals);

   return SCIP_OKAY;
}

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

/* get the "best" master variable among a set of candidates w.r.t. pseudo costs and remove it from the candidate list */
static
SCIP_RETCODE getAndRemoveBestMastercand(
      SCIP*                   scip,
      SCIP_VAR**              mastercands,
      int*                    nmastercands,
      SCIP_Real*              activities,
      int*                    blocknr,
      SCIP_VAR**              mastervar,
      int*                    violchange
      )
{
   SCIP* origprob;
   SCIP_VARDATA* vardata;
   int index;

   int i;
   int tmpviolchange;

   /* get original problem */
   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   index = -1;
   *mastervar = NULL;
   *violchange = SCIPinfinity(scip);
   for( i = *nmastercands - 1; i >= 0; i-- )
   {
      vardata = SCIPvarGetData(mastercands[i]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_MASTER);
      assert(!vardata->data.mastervardata.isray); /* TODO: handle rays */
      assert(!vardata->blocknr >= 0); /* TODO: handle copied original variables and linking variables */

//      SCIPdebugMessage("mastercand %d, pricingprob %d, blocknr %d, nidentblocks %d\n",
//            i, vardata->blocknr, blocknr[vardata->blocknr], GCGrelaxGetNIdenticalBlocks(origprob, vardata->blocknr));

      /* ignore the master variable if the corresponding block is already full */
      if( blocknr[vardata->blocknr] < GCGrelaxGetNIdenticalBlocks(origprob, vardata->blocknr) )
      {
         tmpviolchange = getViolationChange(scip, activities, mastercands[i]);
         if( tmpviolchange < *violchange || *mastervar == NULL )
         {
            index = i;
            *mastervar = mastercands[i];
            *violchange = tmpviolchange;
         }
      }
   }

   assert(index != -1);
   assert(*mastervar != NULL);
   assert(*violchange != SCIPinfinity(scip));

   /* remove selected variable from the candidate list */
   for( i = index; i < *nmastercands - 1; i++ )
   {
      mastercands[i] = mastercands[i + 1];
   }
   (*nmastercands)--;

   return SCIP_OKAY;
}

/* remove master candidates whose blocks are already full */
static
SCIP_RETCODE cleanMastercands(
      SCIP*                   origprob,
      SCIP_VAR**              mastercands,
      int*                    nmastercands,
      int*                    blocknr
      )
{
   SCIP_VARDATA* vardata;
   
   int i;
   int j;
   
   /* clean the candidate list from master variables whose blocks are full */
   for( j = 0; j < *nmastercands; j++ )
   {
      vardata = SCIPvarGetData(mastercands[j]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_MASTER);

      if( blocknr[vardata->blocknr] >= GCGrelaxGetNIdenticalBlocks(origprob, vardata->blocknr) )
         break;
   }
   for( i = j + 1; i < *nmastercands; i++ )
   {
      vardata = SCIPvarGetData(mastercands[i]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_MASTER);

      if( blocknr[vardata->blocknr] < GCGrelaxGetNIdenticalBlocks(origprob, vardata->blocknr) )
      {
         mastercands[j] = mastercands[i];
         j++;
      }
   }

   *nmastercands = j;
   
   return SCIP_OKAY;
}

/* get the "best" master variable w.r.t. pseudo costs */
static
SCIP_RETCODE getBestMastervar(
      SCIP*                   scip,
      SCIP_Real*              activities,
      int*                    blocknr,
      SCIP_VAR**              mastervar,
      int*                    violchange
      )
{
   SCIP* origprob;
   SCIP_VAR** mastervars;
   SCIP_VARDATA* vardata;
   int nmastervars;

   int i;
//   int j;
   int tmpviolchange;

   /* get original problem */
   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   /* get variable data of the master problem */
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(nmastervars >= 0);

   *mastervar = NULL;
   *violchange = SCIPinfinity(scip);

//   j = nmastervars - 1;
//   do
//   {
//      *mastervar = mastervars[j];
//      *violchange = getViolationChange(scip, activities, *mastervar);
//
//      vardata = SCIPvarGetData(*mastervar);
//      assert(vardata != NULL);
//      assert(vardata->vartype == GCG_VARTYPE_MASTER);
//
//      j--;
//   }
//   while( blocknr[vardata->blocknr] >= GCGrelaxGetNIdenticalBlocks(origprob, vardata->blocknr) );

   for( i = nmastervars - 1; i >= 0; i-- )
   {
      vardata = SCIPvarGetData(mastervars[i]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_MASTER);

      /* TODO: handle copied original variables and linking variables */
      if( vardata->blocknr < 0 )
         continue;

      /* ignore the master variable if the corresponding block is already full */
      if( blocknr[vardata->blocknr] < GCGrelaxGetNIdenticalBlocks(origprob, vardata->blocknr)
            && !vardata->data.mastervardata.isray ) /* TODO: handle rays */
      {
         tmpviolchange = getViolationChange(scip, activities, mastervars[i]);
         if( tmpviolchange < *violchange )
         {
            *mastervar = mastervars[i];
            *violchange = tmpviolchange;
         }
      }
   }

   return SCIP_OKAY;
}

/* update working original solution */
static
SCIP_RETCODE updateOrigsol(
      SCIP*                   origprob,
      SCIP_HEUR*              heur,
      SCIP_SOL*               origsol,
      SCIP_VAR*               mastervar,
      int                     violchange,
      int*                    nviolrows,
      SCIP_Real*              activities,
      int*                    blocknr,
      SCIP_Bool*              allblocksfull,
      SCIP_Bool*              success
      )
{
   SCIP_VARDATA* vardata;
   SCIP_VARDATA* vardata2;
   int npricingprobs;

   int i;

   /* get number of pricing problems */
   npricingprobs = GCGrelaxGetNPricingprobs(origprob);

   /* get master variable data */
   vardata = SCIPvarGetData(mastervar);
   assert(vardata != NULL);
   assert(vardata->vartype == GCG_VARTYPE_MASTER);
   assert(vardata->data.mastervardata.norigvars >= 0);
   assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
   assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);
   assert(!vardata->data.mastervardata.isray);

   /* increase master value by one, i.e. increase solution values in current original solution accordingly */
   if( vardata->blocknr == -1 )
   {
      assert(vardata->data.mastervardata.norigvars == 1);
      assert(vardata->data.mastervardata.origvals[0] == 1.0);

      /* increase the corresponding value */
      SCIP_CALL( SCIPincSolVal(origprob, origsol, vardata->data.mastervardata.origvars[0], vardata->data.mastervardata.origvals[0]) );

      /* try to add original solution to solution pool */
      SCIP_CALL( SCIPtrySol(origprob, origsol, FALSE, TRUE, TRUE, TRUE, success) );
   }
   else
   {
      /* loop over all original variables contained in the current master variable */
      for( i = 0; i < vardata->data.mastervardata.norigvars; i++ )
      {
         assert(!SCIPisZero(origprob, vardata->data.mastervardata.origvals[i]));

         /* get the right original variable */
         vardata2 = SCIPvarGetData(vardata->data.mastervardata.origvars[i]);
         assert(vardata2 != NULL);
         assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);

         if(vardata2->blocknr == -2)
            continue;

         assert(vardata2->data.origvardata.pricingvar != NULL);
         vardata2 = SCIPvarGetData(vardata2->data.origvardata.pricingvar);
         assert(vardata2 != NULL);
         assert(vardata2->vartype == GCG_VARTYPE_PRICING);

         /* increase the corresponding value */
         SCIP_CALL( SCIPincSolVal(origprob, origsol, vardata2->data.pricingvardata.origvars[blocknr[vardata->blocknr]], vardata->data.mastervardata.origvals[i]) );
      }

      blocknr[vardata->blocknr]++;

      /* try to add original solution to solution pool */
      SCIP_CALL( SCIPtrySol(origprob, origsol, FALSE, TRUE, TRUE, TRUE, success) );
   }

   /* update number of violated rows and activities array */
   *nviolrows += violchange;
   SCIP_CALL( updateActivities(origprob, activities, mastervar) );

   /* check if all blocks are full */
   *allblocksfull = TRUE;
   for( i = 0; i < npricingprobs; i++ )
   {
      *allblocksfull &= blocknr[i] >= GCGrelaxGetNIdenticalBlocks(origprob, i);
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of primal heuristic
 */


/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#define heurCopyRelaxcolsel NULL

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRelaxcolsel)
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
SCIP_DECL_HEURINIT(heurInitRelaxcolsel)
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
#define heurExitRelaxcolsel NULL


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolRelaxcolsel NULL


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolRelaxcolsel NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRelaxcolsel)
{  /*lint --e{715}*/
   SCIP* origprob;                           /* SCIP structure of original problem  */
   SCIP_HEURDATA* heurdata;                  /* heuristic's data                    */
   SCIP_SOL* origsol;                        /* working original solution           */
   SCIP_VAR* mastervar;
   SCIP_VAR** mastercands;                   /* master variables which are considered first for increasing           */
   SCIP_Real* activities;                    /* for each master LP row, activity of current master solution          */
   int* blocknr;                             /* for each pricing problem, block we are currently working in          */
   SCIP_Bool allblocksfull;                  /* indicates if all blocks are full, i.e. all convexity constraints are satisfied */
   SCIP_Bool success;
   int minnewcols;                           /* minimum number of new columns necessary for calling the heuristic    */
   int nlprows;
   int nmastercands;
   int nmastervars;
   int npricingprobs;
   int nviolrows;
   int violchange;

   int i;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );
   assert(SCIPhasCurrentNodeLP(scip));

   /* get original problem */
   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal relaxation solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* get variable data of the master problem */
   SCIP_CALL( SCIPgetVarsData(scip, NULL, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(nmastervars >= 0);

   /* calculate minimum number of new columns necessary for calling the heuristic;
    * this number is influenced by how successful the heuristic was in the past */
   minnewcols = heurdata->mincolumns * (int) 1.0 * ((1.0 + SCIPheurGetNCalls(heur)) / (1.0 + SCIPheurGetNBestSolsFound(heur)));

   /* if there are not enough new columns since last call, abort heuristic */
   if( nmastervars - heurdata->lastncols < minnewcols )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("Executing GCG relaxation based column selection heuristic (nmastervars = %d) ...\n", nmastervars);

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
   nlprows = SCIPgetNLPRows(scip);
   assert( nlprows >= 0);

   /* get memory for working original solution and row activities */
   SCIP_CALL( SCIPcreateSol(origprob, &origsol, heur) );
   SCIP_CALL( SCIPallocBufferArray(scip, &activities, nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mastercands, 1) );
   nmastercands = 0;

   success = FALSE;

   /* initialize working original solution as transformation of rounded down master LP solution
    * and get the candidate master variables for rounding up */
   SCIP_CALL( initializeOrigsol(scip, heur, origsol, &mastercands, &nmastercands, activities, blocknr, &success) );

   /* first, loop over all candidates for rounding up */
   while( nmastercands > 0 && !allblocksfull && !success )
   {
      SCIP_CALL( getAndRemoveBestMastercand(scip, mastercands, &nmastercands, activities, blocknr, &mastervar, &violchange) );
      SCIP_CALL( updateOrigsol(origprob, heur, origsol, mastervar, violchange, &nviolrows, activities, blocknr, &allblocksfull, &success) );
      SCIP_CALL( cleanMastercands(origprob, mastercands, &nmastercands, blocknr) );

   }

   /* then, consider all master variables for increasing */
   while( !allblocksfull && !success )
   {
      SCIP_CALL( getBestMastervar(scip, activities, blocknr, &mastervar, &violchange) );
      if( mastervar == NULL )
         break;
      SCIP_CALL( updateOrigsol(origprob, heur, origsol, mastervar, violchange, &nviolrows, activities, blocknr, &allblocksfull, &success) );
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
   SCIPfreeBufferArray(scip, &mastercands);

   heurdata->lastncols = nmastervars;

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the relaxation based column selection primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRelaxcolsel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create relaxation based column selection primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyRelaxcolsel, heurFreeRelaxcolsel, heurInitRelaxcolsel, heurExitRelaxcolsel,
         heurInitsolRelaxcolsel, heurExitsolRelaxcolsel, heurExecRelaxcolsel,
         heurdata) );

   /* add relaxation based column selection primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/relaxcolsel/mincolumns",
         "minimum number of columns to regard in the master problem",
         &heurdata->mincolumns, FALSE, DEFAULT_MINCOLUMNS, 1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}