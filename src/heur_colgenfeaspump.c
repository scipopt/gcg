/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_colgenfeaspump.c
 * @ingroup PRIMALHEURISTICS
 * @brief  column generation based feasibility pump primal heuristic
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* toggle debug mode */
//#define SCIP_DEBUG

#include <assert.h>
#include <string.h>

#include "heur_colgenfeaspump.h"
#include "pricer_gcg.h"
#include "pub_gcgvar.h"
#include "relax_gcg.h"
#include "sepa_master.h"

#include "scip/cons_linear.h"
#include "scip/lpi.h"
#include "scip/scipdefplugins.h"


#define HEUR_NAME             "colgenfeaspump"
#define HEUR_DESC             "column generation based feasibility pump"
#define HEUR_DISPCHAR         'G'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXLPITERQUOT    0.01   /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS     1000   /**< additional number of allowed LP iterations */
#define DEFAULT_CYCLELENGTH      20
#define DEFAULT_MAXLOOPS         100    /**< maximal number of pumping rounds */
#define DEFAULT_MAXSTALLLOOPS    10     /**< maximal number of pumping rounds without fractionality improvement (-1: no limit) */
#define DEFAULT_OBJFACTOR        0.95
#define DEFAULT_SHIFTRATE        0.05   /**< percentage of variables to be shifted in case of a 1-cycle */

#define MINLPITER                5000  /**< minimal number of LP iterations allowed in each LP solving call */
#define BIG_M                    100   /**< penalty coefficient for objective funtion */




/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   /* parameters */
   int               cyclelength;
   SCIP_Real         maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   int               maxlpiterofs;       /**< additional number of allowed LP iterations */
   int               maxloops;           /**< maximal number of pumping rounds */
   int               maxstallloops;      /**< maximal number of pumping rounds without fractionality improvement (-1: no limit) */
   SCIP_Real         objfactor;
   SCIP_Real         shiftrate;          /**< percentage of variables to be shifted in case of a 1-cycle */

   /* statistics */
   SCIP_Longint      nlpiterations;      /**< number of LP iterations used in this heuristic */
   int               nsuccess;           /**< number of runs that produced at least one feasible solution */

   /* data */
   int*              masterlocksup;
   int*              masterlocksdown;
   int               nvars;
};




/*
 * Local methods
 */

/*
 * Methods for LP solving, using the LP interface of SCIP
 */

/** copy the master LP to a new SCIP_LPI instance */
static
SCIP_RETCODE initializeLP(
   SCIP*                 scip,
   SCIP_LPI*             divinglp,
   int**                 col2idx,
   int**                 idx2col
   )
{
   SCIP* masterprob;

   SCIP_ROW** masterrows;
   int nmasterrows;
   SCIP_COL** mastercols;
   int nmastercols;

   SCIP_VAR** mastervars;
   int nmastervars;

   /* LP row data */
   SCIP_Real* lhs;
   SCIP_Real* rhs;
   char** rownames;

   /* LP column data */
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   char** colnames;
   int ncolnonz;
   int* colbeg;
   int* colind;
   SCIP_Real* colval;

   int i;
   int j;

   assert(divinglp != NULL);

   /* get master problem */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   /* get master LP rows and columns */
   masterrows = SCIPgetLPRows(masterprob);
   nmasterrows = SCIPgetNLPRows(masterprob);
   assert(nmasterrows >= 0);
   mastercols = SCIPgetLPCols(masterprob);
   nmastercols = SCIPgetNLPCols(masterprob);
   assert(mastercols != NULL);
   assert(nmastercols >= 0);

   /* get master variables' data */
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(nmastercols <= nmastervars);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, col2idx, nmastervars) );
   SCIP_CALL( SCIPallocBufferArray(scip, idx2col, nmastercols) );

   /* allocate memory for rows */
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nmasterrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nmasterrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rownames, nmasterrows) );

   /* get the master LP rows and store them to the new LP
    * (columns are added below) */
   for( i = 0; i < nmasterrows; ++i )
   {
      SCIP_ROW* row;

      /* get row */
      row = masterrows[i];
      assert(row != NULL);

      /* store lhs, rhs and name */
      lhs[i] = SCIProwGetLhs(row);
      rhs[i] = SCIProwGetRhs(row);
      rownames[i] = (char*) SCIProwGetName(row);
   }

   /* store the rows in the new LP */
   SCIP_CALL( SCIPlpiAddRows(divinglp, nmasterrows, lhs, rhs, rownames, 0, NULL, NULL, NULL) );

   /* free memory for rows */
   SCIPfreeBufferArray(scip, &lhs);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &rownames);

   /* allocate memory for columns */
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, nmastercols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nmastercols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nmastercols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colnames, nmastercols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colbeg, nmastercols) );
   colind = NULL;
   colval = NULL;
   ncolnonz = 0;

   /* initialize mapping from variable index to diving lp index */
   for( i = 0; i < nmastervars; ++i)
      (*col2idx)[i] = -1;

   /* copy the master LP columns */
   for( i = 0; i < nmastercols; ++i)
   {
      SCIP_COL* col;
      int colidx;
      SCIP_ROW** colrows;
      SCIP_Real* colvals;
      int ncolrows;

      /* get column and its index in the LP */
      col = mastercols[i];
      assert(col != NULL);
      colidx = SCIPcolGetIndex(col);
      assert(colidx != -1);

      /* get the rows which contain this column */
      ncolrows = SCIPcolGetNNonz(col);
      assert(ncolrows >= 0);
      if( ncolrows > 0 )
      {
         colrows = SCIPcolGetRows(col);
         colvals = SCIPcolGetVals(col);
         assert(colrows != NULL);
         assert(colvals != NULL);
      }

      /* get objective coefficient, lower bound and upper bound, and name */
      obj[i] = SCIPcolGetObj(col);
      lb[i] = SCIPcolGetLb(col);
      ub[i] = SCIPcolGetUb(col);
      colnames[i] = (char*) SCIPvarGetName(SCIPcolGetVar(col));
      colbeg[i] = ncolnonz;

      /* map indices */
      (*col2idx)[colidx] = i;
      (*idx2col)[i] = colidx;

      /* allocate new memory */
      if( ncolrows > 0 )
      {
         if( ncolnonz == 0 )
         {
            assert(colind == NULL);
            assert(colval == NULL);
            SCIP_CALL( SCIPallocBufferArray(scip, &colind, ncolrows) );
            SCIP_CALL( SCIPallocBufferArray(scip, &colval, ncolrows) );
         }
         else
         {
            assert(colind != NULL);
            assert(colval != NULL);
            SCIP_CALL( SCIPreallocBufferArray(scip, &colind, ncolnonz + ncolrows) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &colval, ncolnonz + ncolrows) );
         }
      }

      /* store each entry in the coefficient matrix */
      for( j = 0; j < ncolrows; ++j )
      {
         colind[ncolnonz + j] = SCIProwGetIndex(colrows[j]);
         colval[ncolnonz + j] = colvals[j];
      }
      ncolnonz += ncolrows;
   }

   /* store the columns to the new LP */
   SCIP_CALL( SCIPlpiAddCols(divinglp, nmastercols, obj, lb, ub, colnames, ncolnonz, colbeg, colind, colval) );

   /* free memory for columns */
   SCIPfreeBufferArray(scip, &obj);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &colnames);
   SCIPfreeBufferArray(scip, &colbeg);
   SCIPfreeBufferArrayNull(scip, &colind);
   SCIPfreeBufferArrayNull(scip, &colval);

   return SCIP_OKAY;
}

/** add new variables (columns) to the copied master LP */
static
SCIP_RETCODE addVariables(
   SCIP*                scip,
   SCIP_LPI*            divinglp,
   int**                col2idx,
   int**                idx2col,
   SCIP_VAR**           newvars,
   int                  nnewvars
   )
{
   SCIP* masterprob;

   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_ROW** mastercuts;
   int nmastercuts;

   int nrows;
   int ncols;

   /* LP column data */
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   char** names;
   int nnonz;
   int* beg;
   int* ind;
   SCIP_Real* val;

   int i;
   int j;
   int k;

   /* do not try to add variables if there aren't any */
   if ( nnewvars == 0 )
      return SCIP_OKAY;

   assert(divinglp != NULL);
   assert(col2idx != NULL);
   assert(idx2col != NULL);
   assert(newvars != NULL);
   assert(nnewvars >= 1);

   /* get master problem */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   /* get linear master constraints and cuts */
   masterconss = GCGrelaxGetMasterConss(scip);
   nmasterconss = GCGrelaxGetNMasterConss(scip);
   assert(masterconss != NULL);
   assert(nmasterconss >= 0);
   mastercuts = GCGsepaGetMastercuts(masterprob);
   nmastercuts = GCGsepaGetNMastercuts(masterprob);
   assert(mastercuts != NULL);
   assert(nmastercuts >= 0);

   /* get number of rows in copied LP */
   SCIP_CALL( SCIPlpiGetNRows(divinglp, &nrows) );
   assert(nrows <= nmasterconss + nmastercuts + GCGrelaxGetNPricingprobs(scip));

   /* get number of columns in copied LP */
   SCIP_CALL( SCIPlpiGetNCols(divinglp, &ncols) );
   assert(ncols >= 0);

   /* reallocate memory for mappings */
   SCIP_CALL( SCIPreallocBufferArray(scip, col2idx, SCIPgetNVars(masterprob)) );
   SCIP_CALL( SCIPreallocBufferArray(scip, idx2col, ncols + nnewvars) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, nnewvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nnewvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nnewvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &names, nnewvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &beg, nnewvars) );
   ind = NULL;
   val = NULL;
   nnonz = 0;

   /* for each new master variable, get the coefficients in the master constraints and master cuts
    * and add the to the diving LP */
   for( i = 0; i < nnewvars; ++i )
   {
      SCIP_VAR* newvar;
      int varidx;
      int block;

      newvar = newvars[i];
      assert(newvar != NULL);
      assert(GCGvarIsMaster(newvar));

      varidx = SCIPvarGetProbindex(newvars[i]);
      assert(varidx != -1);

      block = GCGvarGetBlock(newvar);
      assert(block >= 0 && block < GCGrelaxGetNPricingprobs(scip));

      /* get objective coefficient, lower bound and upper bound, and name */
      obj[i] = SCIPvarGetObj(newvar);
      lb[i] = SCIPvarGetLbLocal(newvar);
      ub[i] = SCIPvarGetUbLocal(newvar);
      names[i] = (char*) SCIPvarGetName(newvar);
      beg[i] = nnonz;

      /* allocate new memory */
      if( nnonz == 0 )
      {
         assert(ind == NULL);
         assert(val == NULL);
         SCIP_CALL( SCIPallocBufferArray(scip, &ind, nrows) );
         SCIP_CALL( SCIPallocBufferArray(scip, &val, nrows) );
      }
      else
      {
         assert(ind != NULL);
         assert(val != NULL);
         SCIP_CALL( SCIPreallocBufferArray(scip, &ind, nnonz + nrows) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &val, nnonz + nrows) );
      }

      /* get coefficients for master constraints */
      for( j = 0; j < nmasterconss; ++j )
      {
         SCIP_CONS* cons;
         SCIP_VAR** consvars;
         SCIP_Real* consvals;
         int nconsvars;

         SCIP_ROW* row;
         int idx;

         /* get constraint;
          * @todo: what if the constraint has been upgraded? */
         cons = masterconss[j];
         assert(SCIPconsGetHdlr(cons) != NULL);
         assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "linear") == 0);

         /* get entries in the constraint */
         nconsvars = SCIPgetNVarsLinear(scip, cons);
         assert(nconsvars >= 0);
         if( nconsvars > 0 )
         {
            consvars = SCIPgetVarsLinear(scip, cons);
            consvals = SCIPgetValsLinear(scip, cons);
            assert(consvars != NULL);
            assert(consvals != NULL);
         }

         /* search variable in the constraint */
         for( k = nconsvars - 1; k >= 0; --k )
         {
            if( consvars[k] == newvars[i] )
               break;
         }

         /* if the new variable is in the constraint, add a coefficient for the LP */
         if( k != -1 && !SCIPisZero(scip, consvals[k]) )
         {
            row = SCIPgetRowLinear(scip, cons);
            assert(row != NULL);
            idx = SCIProwGetLPPos(row);

            ind[nnonz] = idx;
            val[nnonz] = consvals[k];
            ++nnonz;
         }
      }

      /* get coefficient in the right convexity constraint */
      if( !GCGmasterVarIsRay(newvar) )
      {
         SCIP_CONS* cons;
         SCIP_ROW* row;
         int idx;

         cons = GCGrelaxGetConvCons(scip, block);
         assert(cons != NULL);
         row = SCIPgetRowLinear(scip, cons);
         assert(row != NULL);
         idx = SCIProwGetLPPos(row);
         assert(idx >= 0);

         ind[nnonz] = idx;
         val[nnonz] = 1.0;
         ++nnonz;
      }

      /* get coefficients for the master cuts */
      for( j = 0; j < nmastercuts; ++j )
      {
         SCIP_ROW* row;
         SCIP_COL** rowcols;
         SCIP_Real* rowvals;
         int nrowcols;

         SCIP_VAR* var;
         int idx;

         /* get entries in the cut */
         row = mastercuts[j];
         rowcols = SCIProwGetCols(row);
         assert(nrowcols >= 0);
         if( nrowcols > 0 )
         {
            rowvals = SCIProwGetVals(row);
            nrowcols = SCIProwGetNNonz(row);
            assert(rowcols != NULL);
            assert(rowvals != NULL);
         }

         /* search variable in the cut */
         for( k = nrowcols - 1; k >= 0; --k)
         {
            var = SCIPcolGetVar(rowcols[k]);
            assert(var != NULL);
            if( var == newvars[i] )
               break;
         }

         /* if the new variable is in the cut, add a coefficient for the LP */
         if( k != -1 && !SCIPisZero(scip, rowvals[k]) )
         {
            idx = SCIProwGetLPPos(row);

            ind[nnonz] = idx;
            val[nnonz] = rowvals[k];
            ++nnonz;
         }
      }

      /* map variable index to column index in diving LP */
      (*col2idx)[varidx] = ncols + i;
      (*idx2col)[ncols + i] = varidx;
   }

   /* add new columns to the diving LP */
   SCIP_CALL( SCIPlpiAddCols(divinglp, nnewvars, obj, lb, ub, names, nnonz, beg, ind, val) );

   /* solve the LP again */
   SCIP_CALL( SCIPlpiSolvePrimal(divinglp) );

   /* free memory */
   SCIPfreeBufferArray(scip, &obj);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &names);
   SCIPfreeBufferArray(scip, &beg);
   SCIPfreeBufferArrayNull(scip, &ind);
   SCIPfreeBufferArrayNull(scip, &val);

   return SCIP_OKAY;
}

/** set new objective coefficients for the LP columns */
static
SCIP_RETCODE setObjectives(
   SCIP*                scip,
   SCIP_LPI*            divinglp,
   int*                 idx2col,
   SCIP_Real*           objectives
   )
{
   SCIP* masterprob;
   SCIP_VAR** mastervars;
   int nmastervars;

   int ncols;
   int* ind;
   SCIP_Real* obj;

   int i;

   assert(divinglp != NULL);
   assert(idx2col != NULL);
   assert(objectives != NULL);

   /* get master problem */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   /* get master variables */
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   /* get number of LP columns */
   SCIP_CALL ( SCIPlpiGetNCols(divinglp, &ncols) );
   assert(ncols <= nmastervars);

   if( ncols == 0 )
      return SCIP_OKAY;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &ind, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, ncols) );

   /* for each master variable which is in the LP, get the new objective */
   for( i = 0; i < ncols; ++i )
   {
      int idx;

      idx = idx2col[i];
      assert(idx >= 0 && idx < nmastervars);

      ind[i] = i;
      obj[i] = objectives[idx];
   }

   /* change objectives */
   SCIP_CALL( SCIPlpiChgObj(divinglp, ncols, ind, obj) );

   /* free memory */
   SCIPfreeBufferArray(scip, &ind);
   SCIPfreeBufferArray(scip, &obj);

   return SCIP_OKAY;
}

/** solve the LP, and store the result into a SCIP_SOL */
static
SCIP_RETCODE solveLP(
   SCIP*                scip,
   SCIP_LPI*            divinglp,
   int*                 col2idx,
   SCIP_HEUR*           heur,
   SCIP_SOL**           lpsol,
   SCIP_Bool*           solved
   )
{
   SCIP* masterprob;

   SCIP_VAR** mastervars;
   SCIP_Real* solvals;
   int nmastervars;

   SCIP_Real* primsol;
   SCIP_Real objval;
   int ncols;

   int i;

   assert(divinglp != NULL);
   assert(col2idx != NULL);

   /* get master problem */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   /* get master variables' data */
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPlpiGetNCols(divinglp, &ncols) );
   assert(ncols <= nmastervars);

   /* free previous lp solution */
   if( *lpsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(masterprob, lpsol));
      *lpsol = NULL;
   }

   /* solve the LP */
   SCIP_CALL( SCIPlpiSolvePrimal(divinglp) );

   /* if the LP was solved, store the obtained solution */
   if( *solved
         && !SCIPlpiIsPrimalInfeasible(divinglp)
         && !SCIPlpiIsPrimalUnbounded(divinglp) )
   {
      SCIP_CALL( SCIPcreateSol(masterprob, lpsol, heur) );

      /* allocate memory for storing the solution values */
      SCIP_CALL( SCIPallocBufferArray(scip, &primsol, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nmastervars) );

      /* get solution and objective value */
      SCIPlpiGetSol(divinglp, &objval, primsol, NULL, NULL, NULL);

      /* get solution values for the master variables */
      for( i = 0; i < nmastervars; ++i)
      {
         int idx;

         idx = col2idx[i];
         assert(-1 <= idx && idx < ncols);
         if( idx != -1 )
            solvals[i] = primsol[idx];
         else
            solvals[i] = 0.0;
      }

      /* store the solution values */
      SCIP_CALL( SCIPsetSolVals(masterprob, *lpsol, nmastervars, mastervars, solvals) );

      /* free memory */
      SCIPfreeBufferArray(scip, &primsol);
      SCIPfreeBufferArray(scip, &solvals);
   }

   return SCIP_OKAY;
}


/*
 * further local methods
 */


/** for a given solution, calculate the number of fractional variables that should be integral */
static
SCIP_RETCODE getNSolFracs(
   SCIP*                scip,
   SCIP_SOL*            relaxsol,
   int*                 nfracs
   )
{
   SCIP_VAR** vars;
   SCIP_Real* solvals;
   int nvars;

   int i;

   *nfracs = 0;

   /* get variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* get solution values */
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, relaxsol, nvars, vars, solvals) );

   for( i = 0; i < nvars; ++i )
   {
      SCIP_VARTYPE vartype;
      SCIP_Real frac;

      vartype = SCIPvarGetType(vars[i]);
      frac = SCIPfeasFrac(scip, solvals[i]);

      if( vartype <= SCIP_VARTYPE_INTEGER && !SCIPisFeasZero(scip, frac) )
         ++(*nfracs);
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &solvals);

   return SCIP_OKAY;
}


/** solve the subproblems with a distance objective function */
static
SCIP_RETCODE solvePricingProblems(
      SCIP*             scip,
      SCIP_HEURDATA*    heurdata,
      SCIP_Real         alpha,
      SCIP_SOL*         relaxsol,
      SCIP_SOL*         sol,
      SCIP_Real*        pricingobjs
      )
{
   SCIP* masterprob;

   SCIP_VAR** origvars;
   int norigvars;
   SCIP_VAR* origvar;
   SCIP_Real solval;
   SCIP_Real frac;
   SCIP_Real newobjcoeff;
   int idx;
   int nlocksup;
   int nlocksdown;

   int npricingprobs;
   SCIP* pricingprob;
   SCIP_VAR** subvars;
   int nsubvars;
   int nbinvars;
   int nintvars;
   SCIP_SOL* subsol;

   int i;
   int j;
   int v;

   /* get master problem and number of pricing problems */
   masterprob = GCGrelaxGetMasterprob(scip);
   npricingprobs = GCGrelaxGetNPricingprobs(scip);

   /* for each pricing problem, change the objective coefficients and solve it */
   for( i = 0; i < npricingprobs; i++ )
   {
      /* get the pricing problem */
      pricingprob = GCGrelaxGetPricingprob(scip, i);
      assert(pricingprob == NULL || GCGrelaxGetNIdenticalBlocks(scip, i) > 0);
      assert(pricingprob != NULL || GCGrelaxGetNIdenticalBlocks(scip, i) == 0);

      /* Due to identical blocks, it may be that the pricing problem of the current block
       * is represented by another one */
      if( pricingprob != NULL )
      {
         int nidenticalblocks;

         /* get the pricing variables and the number of pricing problems represented by the current problem */
         SCIP_CALL( SCIPgetVarsData(pricingprob, &subvars, &nsubvars, &nbinvars, &nintvars, NULL, NULL) );
         nidenticalblocks = GCGrelaxGetNIdenticalBlocks(scip, i);

         /* The pricing problem may represent a number of other pricing problems
          * (in case of identical blocks); in that case, it has to be solved
          * once for each block */
         for( j = 0; j < nidenticalblocks; j++ )
         {
            /* change objective function values;
             * first, look at the binary and integer variables */
            for( v = 0; v < nbinvars + nintvars; v++ )
            {
               assert(GCGvarIsPricing(subvars[v]));
               origvars = GCGpricingVarGetOrigvars(subvars[v]);
               norigvars = GCGpricingVarGetNOrigvars(subvars[v]);
               assert(j < norigvars);

               /* get corresponding variable in the original problem, its index,
                * relaxation solution value and its fractionality */
               origvar = origvars[j];
               idx = SCIPvarGetProbindex(origvar);
               solval = SCIPgetSolVal(scip, relaxsol, origvar);
               frac = SCIPfeasFrac(scip, solval);

               /* compute the objective coefficient;
                * variables which are already integral, are treated separately */
               if( SCIPisFeasZero(scip, frac) )
               {
                  SCIP_Real lb;
                  SCIP_Real ub;

                  /* variables at their bounds should be kept there */
                  lb = SCIPvarGetLbLocal(origvar);
                  ub = SCIPvarGetUbLocal(origvar);
                  if( SCIPisFeasEQ(scip, solval, lb) )
                     newobjcoeff = BIG_M;
                  else if( SCIPisFeasEQ(scip, solval, ub) )
                     newobjcoeff = -BIG_M;
                  else
                     newobjcoeff = 0.0;
               }
               else
               {
                  /* decide by the number of locks (w.r.t. the master constraints)
                   * in which direction the variable should preferably go */
                  nlocksup = heurdata->masterlocksup[idx];
                  nlocksdown = heurdata->masterlocksdown[idx];

                  if( nlocksup > nlocksdown )
                     newobjcoeff = (SCIP_Real) nlocksup / (SCIP_Real) (nlocksup + nlocksdown)
                                    * (SCIPfeasCeil(scip, solval) - solval);
                  else if ( nlocksdown > nlocksup )
                     newobjcoeff = - (SCIP_Real) nlocksdown / (SCIP_Real) (nlocksup + nlocksdown)
                                    * (solval - SCIPfeasFloor(scip, solval));
                  else
                     newobjcoeff = 0.0;
               }

               /* change the objective coefficient */
               SCIP_CALL( SCIPchgVarObj(pricingprob, subvars[v], newobjcoeff) );
               pricingobjs[idx] = newobjcoeff;

               /* reset the solution value to zero */
               SCIP_CALL( SCIPsetSolVal(scip, sol, origvar, 0.0) );
            }

            /* now, look at continous variables; all of them will get objective coefficient zero */
            for( v = nbinvars + nintvars; v < nsubvars; v++ )
            {
               assert(GCGvarIsPricing(subvars[v]));
               origvars = GCGpricingVarGetOrigvars(subvars[v]);
               norigvars = GCGpricingVarGetNOrigvars(subvars[v]);
               assert(j < norigvars);

               /* get corresponding variable in the original problem and its index */
               origvar = origvars[j];
               idx = SCIPvarGetProbindex(origvar);

               /* change the objective coefficient */
               SCIP_CALL( SCIPchgVarObj(pricingprob, subvars[v], 0.0) );
               pricingobjs[idx] = newobjcoeff;

               /* reset the solution value to zero */
               SCIP_CALL( SCIPsetSolVal(scip, sol, origvar, 0.0) );
            }

            /* solve subproblem for current block */
            SCIP_CALL( SCIPsolve(pricingprob) );
            subsol = SCIPgetBestSol(pricingprob);

            /* set solution values of corresponding block in current working solution */
            for( v = 0; v < nsubvars; v++ )
            {
               assert(GCGvarIsPricing(subvars[v]));
               origvars = GCGpricingVarGetOrigvars(subvars[v]);
               norigvars = GCGpricingVarGetNOrigvars(subvars[v]);
               assert(j < norigvars);

               /* get solution value */
               solval = SCIPgetSolVal(pricingprob, subsol, subvars[v]);

               /* solution values which should be integral may not be integral due to numerics;
                * in that case, round them */
               if( SCIPvarGetType(subvars[v]) != SCIP_VARTYPE_CONTINUOUS )
               {
                  assert(SCIPisEQ(scip, solval, SCIPfloor(scip, solval)));
                  solval = SCIPfloor(scip, solval);
               }

               SCIP_CALL( SCIPsetSolVal(scip, sol, origvars[j], solval) );
            }

            /* free pricing problem s. t. it can be solved again */
            SCIP_CALL( SCIPfreeTransform(pricingprob) );
         }
      }
   }

   return SCIP_OKAY;
}

/** check if there are cycles, i. e. if a solution has already been visited before */
static
SCIP_RETCODE checkCycles(
      SCIP*             scip,
      int               cyclelength,
      int               nloops,
      SCIP_SOL*         sol,
      SCIP_Real         alpha,
      SCIP_SOL**        lastsols,
      SCIP_Real*        lastalphas,
      int*              cycle
      )
{
   SCIP_VAR** vars;
   int nvars;

   SCIP_Real solval1;
   SCIP_Real solval2;

   int i;
   int j;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   // TODO: also take alphas into account

   *cycle = -1;
   for( i = 0; i < MIN(cyclelength, nloops - 1); i++ )
   {
      for( j = 0; j < nvars; j++ )
      {
         solval1 = SCIPgetSolVal(scip, sol, vars[j]);
         solval2 = SCIPgetSolVal(scip, lastsols[i], vars[j]);
         if( !SCIPisFeasEQ(scip, solval1, solval2) )
            break;
      }
      if( j == nvars )
      {
         *cycle = i;
         break;
      }
   }

   return SCIP_OKAY;
}

/** shift a solution in case of a 1-cycle */
// TODO: how to calculate scores; in particular, a reasonable weighting between nLocks and nVarshifts
static
SCIP_RETCODE shiftSol(
      SCIP*             scip,
      SCIP_SOL*         sol,
      SCIP_Real         shiftrate,
      SCIP_Real*        pricingobjs,
      int*              nshifts
      )
{
   SCIP_VAR** vars;
   int nbinvars;
   int nintvars;
   int nvars;

   int maxshifts;
   int* varshifts;
   int score;
   int minscore;
   SCIP_VAR* var;
   SCIP_VAR* pricingvar;
   SCIP_VAR* shiftvar;
   int idx;
   SCIP_Bool increase;

   int i;
   int j;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   SCIP_CALL( SCIPallocBufferArray(scip, &varshifts, nvars) );
   for( i = 0; i < nvars; i++ )
      varshifts[i] = 0;

   *nshifts = 0;

   maxshifts = shiftrate * nvars;
   minscore = INT_MAX;
   shiftvar = NULL;
   for( i = 0; i < maxshifts; i++ )
   {
      for( j = 0; j < nvars; j++ )
      {
         if( pricingobjs[j] == 0 || varshifts[j] == INT_MAX)
            continue;

         var = vars[j];
         assert(GCGvarIsOriginal(var));
         pricingvar = GCGoriginalVarGetPricingVar(var);

         score = pricingobjs[j] > 0 ? SCIPvarGetNLocksUp(pricingvar) + varshifts[j] : SCIPvarGetNLocksDown(pricingvar) + varshifts[j];
         if( score < minscore )
         {
            minscore = score;
            shiftvar = var;
            idx = j;
            increase = pricingobjs[j] > 0;
         }
      }

      if( shiftvar == NULL )
         return SCIP_OKAY;

      if( increase )
         SCIP_CALL( SCIPincSolVal(scip, sol, shiftvar, 1.0) );
      else
         SCIP_CALL( SCIPincSolVal(scip, sol, shiftvar, -1.0) );
      varshifts[idx]++;

//      SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, TRUE, &feasible) );
//      if( !feasible )
//      {
//         if( increase )
//            SCIP_CALL( SCIPincSolVal(scip, sol, shiftvar, -1.0) );
//         else
//            SCIP_CALL( SCIPincSolVal(scip, sol, shiftvar, 1.0) );
//         varshifts[j] = INT_MAX;
//      }
//      else
         (*nshifts)++;
   }

   SCIPfreeBufferArray(scip, &varshifts);

   return SCIP_OKAY;
}




/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#define heurCopyColgenfeaspump NULL

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeColgenfeaspump)
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
#define heurInitColgenfeaspump NULL


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitColgenfeaspump NULL


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolColgenfeaspump)
{  /*lint --e{715}*/
   SCIP* masterprob;
   SCIP_HEURDATA* heurdata;

   SCIP_VAR** origvars;
   int nvars;

   int i;
   int j;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get master problem */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   /* get original variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &origvars, &nvars, NULL, NULL, NULL, NULL) );

   /* allocate memory, initialize heuristic's data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->masterlocksup, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &heurdata->masterlocksdown, nvars) );

   heurdata->nvars = nvars;

   /* for each variable, calculate the number of locks */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* var;
      SCIP_CONS** linkingconss;
      SCIP_Real* coefs;
      int ncoefs;

      heurdata->masterlocksup[i] = 0;
      heurdata->masterlocksdown[i] = 0;

      var = origvars[i];
      assert(GCGvarIsOriginal(var));

      /* get constraints transferred to the master problem
       * in which the variable is contained */
      ncoefs = GCGoriginalVarGetNCoefs(var);
      assert(ncoefs >= 0);
      if( ncoefs > 0)
      {
         linkingconss = GCGoriginalVarGetLinkingCons(var);
         coefs = GCGoriginalVarGetCoefs(var);
         assert(linkingconss != NULL);
         assert(coefs != NULL);
      }

      /* for each constraint, check whether there is a lock */
      for( j = 0; j < ncoefs; ++j )
      {
         SCIP_CONS* cons;
         SCIP_Real coef;
         SCIP_Real lhs;
         SCIP_Real rhs;

         /* get constraint;
          * @todo: what if the constraint has been upgraded? */
         cons = linkingconss[j];
         assert(cons != NULL);
         assert(SCIPconsGetHdlr(cons) != NULL);
         assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "linear") == 0);

         /* get coefficient of variable, lhs and rhs */
         coef = coefs[j];
         lhs = SCIPgetLhsLinear(scip, cons);
         rhs = SCIPgetRhsLinear(scip, cons);

         /* compute the locks */
         if( SCIPisPositive(scip, coef) )
         {
            if( !SCIPisInfinity(scip, -lhs) )
               ++heurdata->masterlocksdown[i];
            if( !SCIPisInfinity(scip, rhs) )
               ++heurdata->masterlocksup[i];
         }
         if( SCIPisNegative(scip, coef) )
         {
            if( !SCIPisInfinity(scip, -lhs) )
               ++heurdata->masterlocksup[i];
            if( !SCIPisInfinity(scip, rhs) )
               ++heurdata->masterlocksdown[i];
         }
      }

//      SCIPdebugMessage("Variable %s: nlocksup=%d, nlocksdown=%d\n", SCIPvarGetName(origvars[i]),
//                            heurdata->masterlocksup[i], heurdata->masterlocksdown[i]);
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolColgenfeaspump)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free memory */
   SCIPfreeMemoryArray(scip, &heurdata->masterlocksup);
   SCIPfreeMemoryArray(scip, &heurdata->masterlocksdown);

   return SCIP_OKAY;
}


/** calculates an adjusted maximal number of LP iterations */
static
SCIP_Longint adjustedMaxNLPIterations(
   SCIP_Longint          maxnlpiterations,   /**< regular maximal number of LP iterations */
   SCIP_Longint          nsolsfound,         /**< total number of solutions found so far by SCIP */
   int                   nstallloops         /**< current number of stalling rounds */
   )
{
   if( nstallloops <= 1 )
   {
      if( nsolsfound == 0 )
         return 4*maxnlpiterations;
      else
         return 2*maxnlpiterations;
   }
   else
      return maxnlpiterations;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecColgenfeaspump)
{  /*lint --e{715}*/
   SCIP* masterprob;
   SCIP_HEURDATA* heurdata;
   SCIP_LPI* divinglp;           /* diving LP */
//   SCIP_LPSOLSTAT lpsolstat;     /* status of the LP solution */
   SCIP_RETCODE retcode;
   SCIP_SOL** lastsols;
   SCIP_SOL* mastersol;
   SCIP_SOL* relaxsol;
   SCIP_SOL* intsol;
   SCIP_SOL* tmpsol;
   SCIP_VAR** mastervars;
   SCIP_VAR** vars;
   SCIP_Bool solved;
   SCIP_Bool success;
   SCIP_Real alpha;
   SCIP_Real objnorm;         /* Euclidean norm of the objective function, used for scaling */
   SCIP_Real scalingfactor;   /* factor to scale the original objective function with */
   SCIP_Real objfactor;
   SCIP_Real* lastalphas;
   SCIP_Real* masterobjs;
   SCIP_Real* pricingobjs;
   SCIP_Real* solvals;
   SCIP_Longint nlpiterations;    /* number of LP iterations done during one pumping round */
   SCIP_Longint maxnlpiterations; /* maximum number of LP iterations fpr this heuristic */
   SCIP_Longint nsolsfound;       /* number of solutions found by this heuristic */
   SCIP_Longint ncalls;           /* number of calls of this heuristic */
   int* col2idx;                  /* mapping from variable/column index to column index in diving LP */
   int* idx2col;                  /* mapping from diving LP column index to master variable index */
   int bestnfracs;
   int cycle;
   int lastiterations;
   int maxloops;
   int maxstallloops;
   int nfracs;
   int nloops;
   int nmastervars;
   int oldnmastervars;
   int npricingprobs;
//   int nshifts;
   int nstallloops;
   int nvars;

   char probname[SCIP_MAXSTRLEN];

   int i;
   int j;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get master problem */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   /* get number of pricing problems */
   npricingprobs = GCGrelaxGetNPricingprobs(scip);

   /* get original variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* get master variables' data */
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetStage(masterprob) > SCIP_STAGE_SOLVING || SCIPgetLPSolstat(masterprob) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("Not executing CG Feaspump: master LP not solved to optimality.\n");
      return SCIP_OKAY;
   }

   assert(SCIPhasCurrentNodeLP(masterprob));

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(masterprob) == SCIPgetNNodes(masterprob) && SCIPgetDepth(masterprob) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* only call the column generation feasibility pump once at the root */
   if( SCIPgetDepth(scip) == 0 && SCIPheurGetNCalls(heur) > 0 )
      return SCIP_OKAY;

   /* @todo for some reason, the heuristic is sometimes called with an invalid relaxation solution;
    *       in that case, don't execute it */
   if( !SCIPisRelaxSolValid(scip) )
   {
      SCIPdebugMessage("not executing colgen feaspump: invalid relaxation solution (should not happen!)\n");
      return SCIP_OKAY;
   }

   /* calculate the maximal number of LP iterations until heuristic is aborted */
   nlpiterations = SCIPgetNLPIterations(scip) + SCIPgetNLPIterations(masterprob);
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = 10*SCIPheurGetNBestSolsFound(heur) + heurdata->nsuccess;
   maxnlpiterations = (SCIP_Longint)((1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxlpiterquot * nlpiterations);
   maxnlpiterations += heurdata->maxlpiterofs;

   /* don't try to dive, if we took too many LP iterations during diving */
   if( heurdata->nlpiterations >= maxnlpiterations )
      return SCIP_OKAY;

   /* at the first root call, allow more iterations if there is no feasible solution yet */
   if( SCIPheurGetNCalls(heur) == 0 && SCIPgetNSolsFound(scip) == 0 && SCIPgetDepth(scip) == 0 )
      maxnlpiterations += nlpiterations;

   /* allow at least a certain number of LP iterations in this dive */
   maxnlpiterations = MAX(maxnlpiterations, heurdata->nlpiterations + MINLPITER);

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("executing Column Generation Feasibility Pump ...\n");

   /* calculate factor by which alpha is decreased */
   if( heurdata->objfactor == 1.0 )
      objfactor = MIN(1.0 - 0.1 / (SCIP_Real)(1 + SCIPgetNBestSolsFound(scip)), 0.999);
   else
      objfactor = heurdata->objfactor;

   /* calculate maximal number of loops */
   maxloops = (heurdata->maxloops == -1 ? INT_MAX : heurdata->maxloops);
   maxstallloops = (heurdata->maxstallloops == -1 ? INT_MAX : heurdata->maxstallloops);

   /* allocate further memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &masterobjs, nmastervars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pricingobjs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );

   /* allocate memory for cycle handling */
   SCIP_CALL( SCIPallocBufferArray(scip, &lastsols, heurdata->cyclelength) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lastalphas, heurdata->cyclelength) );
   for( i = 0; i < heurdata->cyclelength; i++ )
   {
      SCIP_CALL( SCIPcreateSol(scip, &lastsols[i], heur) );
   }

   /* initialize working solutions */
   mastersol = NULL;
   tmpsol = NULL;
   SCIP_CALL( SCIPcreateSol(scip, &relaxsol, heur) );
   SCIP_CALL( SCIPlinkRelaxSol(scip, relaxsol) );
   SCIP_CALL( SCIPcreateSol(scip, &intsol, heur) );

   /* create a copy of the master LP */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_divingLP", SCIPgetProbName(scip));
   SCIP_CALL( SCIPlpiCreate(&divinglp, probname, SCIPgetObjsense(masterprob)) );
   SCIP_CALL( initializeLP(scip, divinglp, &col2idx, &idx2col) );

   /* in debug mode, check whether the copied master LP yields the same solution */
#ifndef NDEBUG
   SCIP_CALL( solveLP(scip, divinglp, col2idx, heur, &mastersol, &solved) );
   SCIP_CALL( GCGrelaxTransformMastersolToOrigsol(scip, mastersol, &tmpsol) );

//   SCIPdebugMessage("objectives: relaxsol=%g, divinglpsol=%g\n",
//      SCIPsolGetOrigObj(relaxsol), SCIPsolGetOrigObj(tmpsol));

   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real val1;
      SCIP_Real val2;

      val1 = SCIPgetSolVal(scip, relaxsol, vars[i]);
      val2 = SCIPgetSolVal(scip, tmpsol, vars[i]);
      if( !SCIPisEQ(scip, val1, val2) )
      {
         SCIPdebugMessage("WARNING: different values for var %s: relaxsol=%g, divinglpsol=%g\n",
               SCIPvarGetName(vars[i]), val1, val2);
      }
   }

   SCIP_CALL( SCIPfreeSol(scip, &tmpsol) );
   tmpsol = NULL;
#endif

   /* lp was solved to optimality */
   solved = TRUE;
//   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;

   /* change objectives in the master problem */
   /* scale distance function and original objective to the same norm */
   objnorm = SCIPgetObjNorm(scip);
   objnorm = MAX(objnorm, 1.0);
   scalingfactor = SQRT((SCIP_Real)(SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip))) / objnorm;

   nfracs = SCIPgetNExternBranchCands(scip);
   bestnfracs = nfracs;
   lastiterations = 0;
   nloops = 0;
   nstallloops = 0;
   alpha = 1.0;
   cycle = -1;

   while( nfracs > 0
      && heurdata->nlpiterations < adjustedMaxNLPIterations(maxnlpiterations, nsolsfound, nstallloops)
      && nloops < maxloops && nstallloops < maxstallloops
      && !SCIPisStopped(scip) )
   {
      SCIP_Longint nlpiterationsleft;
      int iterlimit;

      nloops++;
      alpha *= objfactor;

      SCIPdebugMessage("CG Feasibility Pump loop %d: %d fractional variables (alpha: %.4f, stall: %d/%d)\n",
         nloops, nfracs, alpha, nstallloops, maxstallloops);

      /* try to round diving lp solution */
      SCIP_CALL( SCIPgetSolVals(scip, relaxsol, nvars, vars, solvals) );
      SCIP_CALL( SCIPsetSolVals(scip, intsol, nvars, vars, solvals) );
      SCIP_CALL( SCIProundSol(scip, intsol, &success) );

      /* if the rounded solution is feasible and better, add it to GCG */
      if( success )
      {
         SCIPdebugMessage(" -> found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, relaxsol));
         //         SCIP_CALL( SCIPtrySol(scip, intsol, FALSE, FALSE, FALSE, FALSE, &success) );
//#ifdef SCIP_DEBUG
//         SCIP_CALL( SCIPtrySol(scip, intsol, TRUE, TRUE, TRUE, TRUE, &success) );
//#else
         SCIP_CALL( SCIPtrySol(scip, intsol, FALSE, TRUE, TRUE, TRUE, &success) );
//#endif
         if( success )
         {
            SCIPdebugMessage(" -> solution was feasible and good enough\n");
            *result = SCIP_FOUNDSOL;
         }
      }

      /* solve all pricing problems and store the result in the current working solution */
      SCIPdebugMessage(" -> solving pricing problem\n");
      SCIP_CALL( solvePricingProblems(scip, heurdata, alpha, relaxsol, intsol, pricingobjs) );
      SCIPdebugMessage(" -> new integer solution, obj=%g\n", SCIPgetSolOrigObj(scip, intsol));

      /* check for cycles */
      SCIP_CALL( checkCycles(scip, heurdata->cyclelength, nloops, intsol, alpha, lastsols, lastalphas, &cycle) );

      /* @todo: What to do when a cycle occurs? */
      if( cycle >= 0 )
      {
         SCIPdebugMessage(" -> cycle of length %d detected\n", cycle + 1);
//         SCIP_CALL( shiftSol(scip, sol, heurdata->shiftrate, pricingobjs, &nshifts) );
//         if( nshifts > 0 )
//         {
//            SCIPdebugMessage(" -> %d shiftings performed\n", nshifts);
//         }
//         else
//         {
//            SCIPdebugMessage(" -> no shifting performed - change alpha\n");
//            alpha /= objfactor * objfactor;
//            alpha = MIN(alpha, 1.0);
//            nstallloops++;
//            continue;
//         }
      }
      else if( cycle > 0)
      {
         SCIPdebugMessage(" -> cycle of length %d detected\n", cycle + 1);
//         alpha /= objfactor * objfactor;
//         alpha = MIN(alpha, 1.0);
//         nstallloops++;
//         continue;
      }

      /* try to add obtained pricing solution to the solution pool; if it is feasible, then stop */
      SCIP_CALL( SCIPtrySol(scip, intsol, FALSE, TRUE, FALSE, TRUE, &success) );
      if( success )
      {
         SCIPdebugMessage(" -> solving pricing problem yielded feasible solution.\n");
         *result = SCIP_FOUNDSOL;
         break;
      }
      else
      {
         SCIPdebugMessage(" -> not feasible for the original problem\n");
      }

      /* add new columns to the master problem and diving lp and update master variables array */
      oldnmastervars = nmastervars;
      SCIP_CALL( GCGpricerTransOrigSolToMasterVars(masterprob, intsol) );
      SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &masterobjs, nmastervars) );
      SCIP_CALL( addVariables(scip, divinglp, &col2idx, &idx2col, &mastervars[oldnmastervars], nmastervars - oldnmastervars) );
      SCIPdebugMessage(" -> added %d new master variables\n", nmastervars - oldnmastervars);


      /* compute objective coefficients in master problem:
       * for each original variable, compute its new objective coefficient (according to the distance function)
       * and add it to all master variables in which it is contained */
      for( i = 0; i < nmastervars; ++i )
         masterobjs[i] = 0.0;
      for( i = 0; i < nvars; ++i )
      {
         SCIP_VAR* var;

         SCIP_VAR** origmastervars;
         SCIP_Real* origmastervals;
         int norigmastervars;

         SCIP_Real intval;
         SCIP_Real relaxval;
         SCIP_Real direction;
         SCIP_Real newobjcoeff;

         /* get original variable, its solution value and fractionality */
         var = vars[i];
         intval = SCIPgetSolVal(scip, intsol, var);
         relaxval = SCIPgetSolVal(scip, relaxsol, var);
         direction = intval - relaxval;

         /* get master variables which contain this variable */
         origmastervars = GCGoriginalVarGetMastervars(vars[i]);
         origmastervals = GCGoriginalVarGetMastervals(vars[i]);
         norigmastervars = GCGoriginalVarGetNMastervars(vars[i]);
         assert(origmastervars != NULL);
         assert(origmastervals != NULL);
         assert(norigmastervars >= 0);

         /* variables which stayed integral are treated separately */
         if( SCIPisFeasZero(scip, direction) )
         {
            SCIP_Real lb;
            SCIP_Real ub;

            /* variables at their bounds should be kept there */
            lb = SCIPvarGetLbLocal(var);
            ub = SCIPvarGetUbLocal(var);
            if( SCIPisFeasEQ(scip, intval, lb) )
               newobjcoeff = BIG_M;
            else if( SCIPisFeasEQ(scip, intval, ub) )
               newobjcoeff = -BIG_M;
            else
               newobjcoeff = 0.0;
         }
         else
         {
            if( SCIPisFeasPositive(scip, direction) )
               newobjcoeff = - 1.0;
            else
               newobjcoeff = 1.0;
         }

         for( j = 0; j < norigmastervars; ++j )
         {
            int idx;

            idx = SCIPvarGetProbindex(origmastervars[j]);
            assert(idx >= 0 && idx < nmastervars);
            masterobjs[idx] += origmastervals[j] * newobjcoeff;
         }
      }

      /* set the new objectives */
      SCIP_CALL( setObjectives(scip, divinglp, idx2col, masterobjs) );

      /* the LP with the new (distance) objective is solved */
      nlpiterations = SCIPgetNLPIterations(scip) + SCIPgetNLPIterations(masterprob) + lastiterations;
      nlpiterationsleft = adjustedMaxNLPIterations(maxnlpiterations, nsolsfound, nstallloops) - heurdata->nlpiterations;
      iterlimit = MAX((int)nlpiterationsleft, MINLPITER);
      SCIPdebugMessage(" -> solve LP with iteration limit %d\n", iterlimit);

      SCIP_CALL( SCIPlpiSetIntpar(divinglp, SCIP_LPPAR_LPITLIM, iterlimit) );
      retcode = solveLP(scip, divinglp, col2idx, heur, &mastersol, &solved);

      /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
       * Hence in optimized mode, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
       */
      if( retcode != SCIP_OKAY )
      {
         SCIPwarningMessage("Error while solving LP in Colgen Feaspump heuristic; LP solve terminated with code <%d>\n", retcode);
         SCIPwarningMessage("This does not affect the remaining solution procedure --> continue\n");
      }

      /* update iteration count;
       * @todo: the method yields an int and not a long int */
      SCIP_CALL( SCIPlpiGetIterations(divinglp, &lastiterations) );
      heurdata->nlpiterations += lastiterations;
      SCIPdebugMessage(" -> number of iterations: %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", solved=%d\n",
            heurdata->nlpiterations, adjustedMaxNLPIterations(maxnlpiterations, nsolsfound, nstallloops), solved);

      /* check whether LP was solved to optimality */
      if( !solved )
         break;

      /* swap the last solutions, i.e. store the pricing solution into the lastsols array */
      tmpsol = lastsols[heurdata->cyclelength-1];
      for( i = heurdata->cyclelength-1; i > 0; i-- )
      {
         lastsols[i] = lastsols[i-1];
         lastalphas[i] = lastalphas[i-1];
      }
      lastsols[0] = intsol;
      lastalphas[0] = alpha;
      intsol = tmpsol;

      /* translate sol into original variable space
       * and check for improvement in number of fractionals */
      SCIP_CALL( SCIPfreeSol(scip, &relaxsol) );
      SCIP_CALL( GCGrelaxTransformMastersolToOrigsol(scip, mastersol, &relaxsol) );
      SCIP_CALL( getNSolFracs(scip, relaxsol, &nfracs) );
      if( nfracs < bestnfracs )
      {
         bestnfracs = nfracs;
         nstallloops = 0;
      }
      else
         nstallloops++;

      SCIPdebugMessage(" -> loop finished: %d fractional variables (stall: %d/%d, iterations: %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT")\n",
            nfracs, nstallloops, maxstallloops, heurdata->nlpiterations, adjustedMaxNLPIterations(maxnlpiterations, nsolsfound, nstallloops));
   }

   /* try final solution, if no more fractional variables are left */
   if( nfracs == 0 && solved )
   {
      success = FALSE;

      SCIPdebugMessage("colgen feaspump found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, relaxsol));

//      SCIP_CALL( SCIPtrySol(scip, relaxsol, FALSE, FALSE, FALSE, FALSE, &success) );
//#ifdef SCIP_DEBUG
//      SCIP_CALL( SCIPtrySol(scip, relaxsol, TRUE, TRUE, TRUE, TRUE, &success) );
//#else
      SCIP_CALL( SCIPtrySol(scip, relaxsol, FALSE, TRUE, TRUE, TRUE, &success) );
//#endif

      if( success )
      {
         SCIPdebugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* free diving LP */
   SCIP_CALL( SCIPlpiFree(&divinglp) );

   /* free memory */
   SCIPfreeBufferArray(scip, &masterobjs);
   SCIPfreeBufferArray(scip, &pricingobjs);
   SCIPfreeBufferArray(scip, &solvals);
   SCIPfreeBufferArray(scip, &col2idx);
   SCIPfreeBufferArray(scip, &idx2col);

   /* free working solutions */
   if( mastersol != NULL )
      SCIP_CALL( SCIPfreeSol(masterprob, &mastersol) );
   if( relaxsol != NULL )
      SCIP_CALL( SCIPfreeSol(scip, &relaxsol) );
   SCIP_CALL( SCIPfreeSol(scip, &intsol) );

   /* free memory for cycle handling */
   for( i = 0; i < heurdata->cyclelength; i++ )
   {
      SCIPfreeSol(scip, &lastsols[i]);
   }
   SCIPfreeBufferArray(scip, &lastsols);
   SCIPfreeBufferArray(scip, &lastalphas);

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the colgenfeaspump primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurColgenfeaspump(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create colgenfeaspump primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyColgenfeaspump,
         heurFreeColgenfeaspump, heurInitColgenfeaspump, heurExitColgenfeaspump,
         heurInitsolColgenfeaspump, heurExitsolColgenfeaspump, heurExecColgenfeaspump,
         heurdata) );

   /* add colgenfeaspump primal heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/"HEUR_NAME"/maxlpiterquot",
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &heurdata->maxlpiterquot, FALSE, DEFAULT_MAXLPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/"HEUR_NAME"/maxlpiterofs",
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, FALSE, DEFAULT_MAXLPITEROFS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/cyclelength",
         "maximum length of cycles to be checked explicitly in each round",
         &heurdata->cyclelength, TRUE, DEFAULT_CYCLELENGTH, 1, 100, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxloops",
         "maximal number of pumping rounds (-1: no limit)",
         &heurdata->maxloops, TRUE, DEFAULT_MAXLOOPS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxstallloops",
         "maximal number of pumping rounds without fractionality improvement (-1: no limit)",
         &heurdata->maxstallloops, TRUE, DEFAULT_MAXSTALLLOOPS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/objfactor",
         "factor by which the regard of the objective is decreased in each round",
         &heurdata->objfactor, FALSE, DEFAULT_OBJFACTOR, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/shiftrate",
         "percentage of variables to be shifted in case of a 1-cycle",
         &heurdata->shiftrate, TRUE, DEFAULT_SHIFTRATE, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
