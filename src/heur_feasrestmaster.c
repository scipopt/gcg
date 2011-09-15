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

/**@file   heur_feasrestmaster.c
 * @ingroup PRIMALHEURISTICS
 * @brief  feasibility restricted master primal heuristic
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* toggle debug mode */
#define SCIP_DEBUG

#include <assert.h>

#include "heur_feasrestmaster.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/pub_misc.h"


#define HEUR_NAME             "feasrestmaster"
#define HEUR_DESC             "heuristic that regards only master variables that cause the least master infeasibility"
#define HEUR_DISPCHAR         'V'
#define HEUR_PRIORITY         100
//#define HEUR_FREQ             10
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          5
#define HEUR_MAXDEPTH         -1
/* TODO: should heuristic be called during the pricing loop or only after solving a node relaxation? */
#define HEUR_TIMING           SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_DURINGPRICINGLOOP
#define HEUR_USESSUBSCIP      TRUE

#define DEFAULT_MAXNODES      5000LL    /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which feasibility restricted master should at least improve the incumbent */
#define DEFAULT_MINNODES      500LL     /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_NODESOFS      500LL     /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_RESTRATE      0.25      /* percentage of master variables that are not fixed to zero           */
#define DEFAULT_USELPROWS     TRUE      /* should subproblem be created out of the rows in the LP rows,
                                         * otherwise, the copy constructor of the constraints handlers are used*/




/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;          /**< maximum number of nodes to regard in the subproblem                 */
   SCIP_Longint          minnodes;          /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Longint          nodesofs;          /**< number of nodes added to the contingent of the total nodes          */
   SCIP_Longint          usednodes;         /**< nodes already used by feasibility restricted master in earlier calls */
   SCIP_Real             restrate;          /**< percentage of master variables that are not fixed to zero           */
   SCIP_Real             minimprove;        /**< factor by which feasibility restricted master should at least improve the incumbent */
   SCIP_Real             nodesquot;         /**< subproblem nodes in relation to nodes of the original problem       */
   SCIP_Bool             uselprows;         /**< should subproblem be created out of the rows in the LP rows?        */
};




/*
 * Local methods
 */

/** sort master variables by their total violation */
static
SCIP_RETCODE sortMastervars(
   SCIP*                 scip,
   SCIP_VAR**            bestmastervars,
   int                   nmastervars
   )
{
   SCIP_COL* col;
   SCIP_ROW** colrows;
   SCIP_Real* colvals;
   SCIP_Real* infeasibilities;
   int ncolrows;

   int v;
   int r;

   SCIP_CALL( SCIPallocBufferArray(scip, &infeasibilities, nmastervars) );

   /* for each master variable, compute the total amount of row infeasibilities it causes */
   for( v = 0; v < nmastervars; v++ )
   {
      infeasibilities[v] = 0;

      col = SCIPvarGetCol(bestmastervars[v]);
      colrows = SCIPcolGetRows(col);
      colvals = SCIPcolGetVals(col);
      ncolrows = SCIPcolGetNLPNonz(col);
      assert(ncolrows == 0 || (colrows != NULL && colvals != NULL));

      for( r = 0; r < ncolrows; r++ )
      {
         SCIP_ROW* row;
         int rowpos;

         row = colrows[r];
         rowpos = SCIProwGetLPPos(row);
         assert(-1 <= rowpos);

         if( rowpos >= 0 && !SCIProwIsLocal(row) )
         {
            assert(SCIProwIsInLP(row));

            infeasibilities[v] += MAX3(0, SCIProwGetLhs(row) - colvals[r], colvals[r] - SCIProwGetRhs(row));
         }
      }
   }

   SCIPsortRealPtr(infeasibilities, (void**) bestmastervars, nmastervars);

   SCIPfreeBufferArray(scip, &infeasibilities);

   return SCIP_OKAY;
}

/** creates a restricted master problem by fixing master variables which are zero */
static
SCIP_RETCODE createSubproblem(
   SCIP*                 scip,               /**< SCIP data structure for master problem                         */
   SCIP*                 restmaster,         /**< SCIP data structure for restricted master problem              */
   SCIP_VAR**            restmastervars,     /**< the variables of the restricted                                */
   SCIP_HASHMAP*         varmapfw,           /**< mapping of master variables to restricted master variables     */
   SCIP_Real             restrate,           /**< percentage of master variables that remain unfixed             */
   SCIP_Bool             uselprows,          /**< should subproblem be created out of the rows in the LP rows?   */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully  */
   )
{
   SCIP_VAR** mastervars;
   SCIP_VAR** bestmastervars;
   SCIP_Bool* fixed;
   int nmastervars;

   int i;

   /* get variable data of the master problem */
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   /* copy master variables array to another array */
   SCIP_CALL( SCIPallocBufferArray(scip, &bestmastervars, nmastervars ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixed, nmastervars ) );
   for( i = 0; i < nmastervars; i++ )
   {
      bestmastervars[i] = mastervars[i];
      fixed[i] = TRUE;
   }

   /* sort master variables by their infeasibilities */
   SCIP_CALL( sortMastervars(scip, bestmastervars, nmastervars) );

   /* mark the master variables with the lowest violation to remain unfixed */
   for( i = 0; i < nmastervars * restrate; i++ )
   {
      int index;

      index = SCIPvarGetProbindex(bestmastervars[i]);
      fixed[index] = FALSE;
   }

   /* create the variables of the restricted master problem */
   for( i = 0; i < nmastervars; i++ )
   {
      /* if master variable is marked, fix it to zero in restricted master */
      if( fixed[i] )
      {
         SCIP_CALL( SCIPcreateVar(restmaster, &restmastervars[i], SCIPvarGetName(mastervars[i]), 0, 0,
               SCIPvarGetObj(mastervars[i]), SCIPvarGetType(mastervars[i]), SCIPvarIsInitial(mastervars[i]),
               SCIPvarIsRemovable(mastervars[i]), NULL, NULL, NULL, NULL, NULL) );
      }

      /* otherwise, just copy the variable to restricted master */
      else
      {
         SCIP_CALL( SCIPcreateVar(restmaster, &restmastervars[i], SCIPvarGetName(mastervars[i]),
               SCIPvarGetLbGlobal(mastervars[i]), SCIPvarGetUbGlobal(mastervars[i]), SCIPvarGetObj(mastervars[i]),
               SCIPvarGetType(mastervars[i]), SCIPvarIsInitial(mastervars[i]), SCIPvarIsRemovable(mastervars[i]),
               NULL, NULL, NULL, NULL, NULL) );
      }

      SCIP_CALL( SCIPaddVar(restmaster, restmastervars[i]) );

      /* insert variable into mapping between master and restricted master */
      SCIP_CALL( SCIPhashmapInsert(varmapfw, mastervars[i], restmastervars[i]) );
   }

   if( uselprows )
   {
      SCIP_ROW** rows;                          /* original scip rows                         */
      int nrows;

      /* get the rows and their number */
      SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

      /* copy all rows to linear constraints */
      for( i = 0; i < nrows; i++ )
      {
         SCIP_CONS* cons;
         SCIP_VAR** consvars;
         SCIP_COL** cols;
         SCIP_Real constant;
         SCIP_Real lhs;
         SCIP_Real rhs;
         SCIP_Real* vals;
         int nnonz;
         int j;

         /* ignore rows that are only locally valid */
         if( SCIProwIsLocal(rows[i]) )
            continue;

         /* get the row's data */
         constant = SCIProwGetConstant(rows[i]);
         lhs = SCIProwGetLhs(rows[i]) - constant;
         rhs = SCIProwGetRhs(rows[i]) - constant;
         vals = SCIProwGetVals(rows[i]);
         nnonz = SCIProwGetNNonz(rows[i]);
         cols = SCIProwGetCols(rows[i]);

         assert( lhs <= rhs );

         /* allocate memory array to be filled with the corresponding subproblem variables */
         SCIP_CALL( SCIPallocBufferArray(restmaster, &consvars, nnonz) );
         for( j = 0; j < nnonz; j++ )
            consvars[j] = restmastervars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

         /* create a new linear constraint and add it to the subproblem */
         SCIP_CALL( SCIPcreateConsLinear(restmaster, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCons(restmaster, cons) );
         SCIP_CALL( SCIPreleaseCons(restmaster, &cons) );

         /* free temporary memory */
         SCIPfreeBufferArray(restmaster, &consvars);
      }
   }

   *success = TRUE;

   SCIPfreeBufferArray(scip, &bestmastervars);
   SCIPfreeBufferArray(scip, &fixed);

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by translating the solution of the restricted master problem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 origprob,           /**< original SCIP data structure                        */
   SCIP*                 scip,               /**< SCIP data structure of master problem               */
   SCIP*                 restmaster,         /**< SCIP structure of restricted master problem         */
   SCIP_VAR**            restmastervars,     /**< the variables of the restricted master problem      */
   SCIP_HEUR*            heur,               /**< RENS heuristic structure                            */
   SCIP_SOL*             restmastersol,      /**< solution of the restricted master problem           */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   SCIP_VAR** mastervars;                    /* the master problem's variables                  */
   int        nvars;
   int        nmastervars;
   SCIP_Real* restmastervals;                /* solution values of the subproblem               */
   SCIP_SOL*  newmastersol;                  /* solution for the master problem                 */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert( origprob != NULL );
   assert( scip != NULL );
   assert( restmaster != NULL );
   assert( restmastervars != NULL );
   assert( restmastersol != NULL );

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(origprob, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert( nmastervars == SCIPgetNOrigVars(restmaster) );

   SCIP_CALL( SCIPallocBufferArray(scip, &restmastervals, nmastervars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(restmaster, restmastersol, nmastervars, restmastervars, restmastervals) );

   /* create new solution for the master problem and translate it to the original problem;
    * TODO: GCG does not recognize that the solution comes from this heuristic */
   SCIP_CALL( SCIPcreateSol(scip, &newmastersol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newmastersol, nmastervars, mastervars, restmastervals) );
   SCIP_CALL( GCGrelaxTransformMastersolToOrigsol(origprob, newmastersol, &newsol) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(origprob, &newsol, FALSE, TRUE, TRUE, TRUE, success) );
   SCIP_CALL( SCIPfreeSol(scip, &newmastersol) );

   SCIPfreeBufferArray(scip, &restmastervals);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#define heurCopyFeasrestmaster NULL

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeFeasrestmaster)
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
SCIP_DECL_HEURINIT(heurInitFeasrestmaster)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* initialize data */
   heurdata->usednodes = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitFeasrestmaster NULL


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolFeasrestmaster NULL


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolFeasrestmaster NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecFeasrestmaster)
{  /*lint --e{715}*/
   SCIP* origprob;                           /* SCIP structure of original problem    */
   SCIP_HEURDATA* heurdata;                  /* heuristic's data                    */
   SCIP_Real timelimit;                      /* timelimit for the subproblem        */
   SCIP_Real memorylimit;
   SCIP_Real cutoff;                         /* objective cutoff for the restricted  master problem        */
   SCIP_Bool discretization;
   SCIP_Bool success;
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the restricted master problem */

   SCIP* restmaster;                         /* SCIP structure of the restricted master problem            */
   SCIP_HASHMAP* varmapfw;                   /* mapping of master variables to restricted master variables */
   SCIP_SOL** restmastersols;
   SCIP_VAR** mastervars;
   SCIP_VAR** restmastervars;
   int nmastervars;
   int nrestmastersols;

   char probname[SCIP_MAXSTRLEN];
   int i;

#ifdef NDEBUG
   SCIP_RETCODE retstat;
#endif

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );
   assert( SCIPhasCurrentNodeLP(scip) );

   /* get original problem */
   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   *result = SCIP_DIDNOTRUN;

   /* this heuristic works only for the discretization approach */
   /* TODO: make heuristic also usable for convexification;
    *       in this case, we need some sort of constraint handler for the restmaster subSCIP */
   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if( !discretization )
      return SCIP_OKAY;

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(origprob));

   /* reward restricted master if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
//      SCIPdebugMessage("skipping GCG feasibility restricted master heuristic: nstallnodes=%"SCIP_LONGINT_FORMAT", minnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(origprob, "limits/time", &timelimit) );
   if( !SCIPisInfinity(origprob, timelimit) )
      timelimit -= SCIPgetSolvingTime(origprob);
   SCIP_CALL( SCIPgetRealParam(origprob, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(origprob, memorylimit) )
      memorylimit -= SCIPgetMemUsed(origprob)/1048576.0;
   if( timelimit < 10.0 || memorylimit <= 0.0 )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPdebugMessage("Executing GCG feasibility restricted master heuristic ...\n");

   *result = SCIP_DIDNOTFIND;

   /* get variable data of the master problem */
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   /* initializing the subproblem */
   SCIP_CALL( SCIPcreate(&restmaster) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(restmaster), SCIPcalcHashtableSize(5 * nmastervars)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &restmastervars, nmastervars) );

   /* include SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(restmaster) );

   /* get name of the master problem and add the string "_restricted" */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_restricted", SCIPgetProbName(scip));

   /* create the subproblem */
   SCIP_CALL( SCIPcreateProb(restmaster, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   //   /* copy all variables */
   //   SCIP_CALL( SCIPcopyVars(masterprob, restmaster, varmapfw, NULL, TRUE) );

   success = FALSE;

   /* create a new problem, which fixes variables with same value in bestsol and LP relaxation */
   SCIP_CALL( createSubproblem(scip, restmaster, restmastervars, varmapfw, heurdata->restrate, heurdata->uselprows, &success) );
   SCIPdebugMessage("feasibility restricted master problem: %d vars, %d cons, success=%u\n", SCIPgetNVars(restmaster), SCIPgetNConss(restmaster), success);

   /* if the lp rows are not used, also copy the constraints */
   if( !heurdata->uselprows )
   {
      SCIP_Bool valid;
      valid = FALSE;

      SCIP_CALL( SCIPcopyConss(scip, restmaster, varmapfw, NULL, TRUE, FALSE, &valid) );
      SCIPdebugMessage("Copying the SCIP constraints was %s complete.\n", valid ? "" : "not ");
   }

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(restmaster, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(restmaster, "display/verblevel", 0) );

   /* set limits for the subproblem */
   SCIP_CALL( SCIPsetLongintParam(restmaster, "limits/stallnodes", nstallnodes) );
   SCIP_CALL( SCIPsetLongintParam(restmaster, "limits/nodes", heurdata->maxnodes) );
   SCIP_CALL( SCIPsetRealParam(restmaster, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(restmaster, "limits/memory", memorylimit) );

   /* forbid recursive call of heuristics solving subMIPs */
   SCIP_CALL( SCIPsetSubscipsOff(restmaster, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(restmaster, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(restmaster, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(scip, "estimate") != NULL )
   {
      SCIP_CALL( SCIPsetIntParam(restmaster, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(scip, "inference") != NULL )
   {
      SCIP_CALL( SCIPsetIntParam(restmaster, "branching/inference/priority", INT_MAX/4) );
   }

   /* disable conflict analysis */
   SCIP_CALL( SCIPsetBoolParam(restmaster, "conflict/useprop", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(restmaster, "conflict/useinflp", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(restmaster, "conflict/useboundlp", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(restmaster, "conflict/usesb", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(restmaster, "conflict/usepseudo", FALSE) );

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* if the subproblem could not be created, free memory and return */
   if( !success )
   {
      *result = SCIP_DIDNOTRUN;
      SCIP_CALL( SCIPfreeTransform(restmaster) );
      for( i = 0; i < nmastervars; i++ )
      {
         SCIP_CALL( SCIPreleaseVar(restmaster, &restmastervars[i]) );
      }
      SCIPfreeBufferArray(scip, &restmastervars);
      SCIP_CALL( SCIPfree(&restmaster) );
      return SCIP_OKAY;
   }

   /* if there is already a solution, add an objective cutoff */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real upperbound;
      assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

      if( !SCIPisInfinity(scip,-1.0*SCIPgetLowerbound(scip)) )
      {
         cutoff = (1-heurdata->minimprove)*SCIPgetUpperbound(scip) + heurdata->minimprove*SCIPgetLowerbound(scip);
      }
      else
      {
         if ( SCIPgetUpperbound ( scip ) >= 0 )
            cutoff = ( 1 - heurdata->minimprove ) * SCIPgetUpperbound ( scip );
         else
            cutoff = ( 1 + heurdata->minimprove ) * SCIPgetUpperbound ( scip );
      }
      cutoff = MIN(upperbound, cutoff);
      SCIP_CALL( SCIPsetObjlimit(restmaster, cutoff) );
   }

   /* solve the restricted master problem */
   /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
    * Hence in optimized mode, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
    */
#ifdef NDEBUG
   retstat = SCIPpresolve(restmaster);
   if( retstat != SCIP_OKAY )
   {
      SCIPwarningMessage("Error while presolving subMIP in GCG restricted master heuristic; restricted master terminated with code <%d>\n",retstat);
   }
#else
   SCIP_CALL( SCIPpresolve(restmaster) );
#endif

   SCIPdebugMessage("presolved restricted master problem: %d vars, %d cons, success=%u\n", SCIPgetNVars(restmaster), SCIPgetNConss(restmaster), success);
   SCIPdebugMessage("solving restricted master problem: nstallnodes=%"SCIP_LONGINT_FORMAT", maxnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, heurdata->maxnodes);

#ifdef NDEBUG
   retstat = SCIPsolve(restmaster);
   if( retstat != SCIP_OKAY )
   {
      SCIPwarningMessage("Error while solving subMIP in GCG feasibility restricted master heuristic; restricted master terminated with code <%d>\n",retstat);
   }
#else
   SCIP_CALL( SCIPsolve(restmaster) );
#endif

   SCIPdebugMessage("GCG feasibility restricted master heuristic: %d feasible solution(s) found.\n", SCIPgetNSols(restmaster));

   /* check, whether a solution was found;
    * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
    */
   nrestmastersols = SCIPgetNSols(restmaster);
   restmastersols = SCIPgetSols(restmaster);
   success = FALSE;
   for( i = 0; i < nrestmastersols && !success; ++i )
   {
      SCIP_CALL( createNewSol(origprob, scip, restmaster, restmastervars, heur, restmastersols[i], &success) );
   }
   if( success )
      *result = SCIP_FOUNDSOL;

   /* free subproblem */
   SCIP_CALL( SCIPfreeTransform(restmaster) );
   for( i = 0; i < nmastervars; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(restmaster, &restmastervars[i]) );
   }
   SCIPfreeBufferArray(scip, &restmastervars);
   SCIP_CALL( SCIPfree(&restmaster) );


   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the feasibility restricted master primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurFeasrestmaster(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create feasibility restricted master primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyFeasrestmaster, heurFreeFeasrestmaster, heurInitFeasrestmaster, heurExitFeasrestmaster,
         heurInitsolFeasrestmaster, heurExitsolFeasrestmaster, heurExecFeasrestmaster,
         heurdata) );

   /* add feasibility restricted master primal heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/feasrestmaster/restrate",
         "percentage of master variables that are not fixed to zero",
         &heurdata->restrate, FALSE, DEFAULT_RESTRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/feasrestmaster/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes,  TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/feasrestmaster/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/feasrestmaster/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/feasrestmaster/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/feasrestmaster/minimprove",
         "factor by which feasibility restricted master should at least improve the incumbent  ",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/feasrestmaster/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   return SCIP_OKAY;
}
