/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_extremepoints.c,v 1.25 2010/01/04 20:35:41 bzfheinz Exp $"

/**@file   heur_extremepoints.c
 * @ingroup PRIMALHEURISTICS
 * @brief  extreme points crossover primal heuristic
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* toggle debug mode */
#define SCIP_DEBUG

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "heur_extremepoints.h"
#include "relax_gcg.h"
#include "gcgplugins.h"
#include "struct_vardata.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"


#define HEUR_NAME             "extremepoints"
#define HEUR_DESC             "heuristic that performs a crossover on the extreme points of a relaxation solution"
#define HEUR_DISPCHAR         'X'
#define HEUR_PRIORITY         -1101500
#define HEUR_FREQ             30
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERPSEUDONODE
#define HEUR_USESSUBSCIP      TRUE

#define DEFAULT_MAXNODES      1000LL        /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINIMPROVE    0.01          /* factor by which crossover should at least improve the incumbent     */
#define DEFAULT_MINNODES      200LL         /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_SIMILARITYRATE 0.5          /* rate by which to decide whether two solutions are "similar"         */
#define DEFAULT_MAXRATEREDNS  3             /* maximum times that groupingrate is allowed to be reduced            */
#define DEFAULT_MINFIXINGRATE 0.4           /* minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_MAXFIXINGRATE 1.000         /* maximum percentage of integer variables that can be fixed           */
#define DEFAULT_ZEROFIXPROB   0.85          /* probability that a zero variable is fixed                           */
#define DEFAULT_NODESOFS      200LL         /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1           /* subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_MAXRUNS       1             /* maximum times the heuristic runs on an LP solution                  */
#define DEFAULT_NUSEDSOLS     2             /* number of solutions that will be taken into account                 */
#define DEFAULT_NWAITINGNODES 200LL         /* number of nodes without incumbent change heuristic should wait      */
#define DEFAULT_RANDOMIZATION TRUE          /* should the choice which sols to take be randomized?                 */
#define DEFAULT_DONTWAITATROOT FALSE        /* should the nwaitingnodes parameter be ignored at the root node?     */
#define DEFAULT_USELPROWS     TRUE           /* should subproblem be created out of the rows in the LP rows,
                                              * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_USEGCG        FALSE         /* should the subproblem be solved with GCG?                           */

#define HASHSIZE_SOLS         11113         /* size of hash table for solution tuples in crossover heuristic       */




/*
 * Data structures
 */
typedef struct SolTuple SOLTUPLE;

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;          /**< maximum number of nodes to regard in the subproblem               */
   SCIP_Longint          minnodes;          /**< minimum number of nodes to regard in the subproblem               */
   SCIP_Longint          nodesofs;          /**< number of nodes added to the contingent of the total nodes        */
   SCIP_Longint          usednodes;         /**< nodes already used by crossover in earlier calls                  */
   SCIP_Real             nodesquot;         /**< subproblem nodes in relation to nodes of the original problem     */

   int                   nusedsols;         /**< number of solutions that will be taken into account               */
   int                   maxruns;           /**< maximum times the heuristic runs on an LP solution                */
   SCIP_Longint          nwaitingnodes;     /**< number of nodes without incumbent change heuristic should wait    */
   unsigned int          nfailures;         /**< number of failures since last successful call                     */
   SCIP_Longint          nextnodenumber;    /**< number of BnB nodes at which crossover should be called next      */
   SCIP_Real             similarityrate;    /**< rate by which to decide whether two solutions are "similar"       */
   int                   maxrateredns;      /**< maximum times that similarityrate is allowed to be reduced        */
   SCIP_Real             minfixingrate;     /**< minimum percentage of integer variables that have to be fixed     */
   SCIP_Real             maxfixingrate;     /**< maximum percentage of integer variables that can be fixed         */
   SCIP_Real             zerofixprob;       /**< probability that a zero variable is fixed                         */
   SCIP_Real             minimprove;        /**< factor by which crossover should at least improve the incumbent   */
   SCIP_Bool             randomization;     /**< should the choice which sols to take be randomized?               */
   SCIP_Bool             dontwaitatroot;    /**< should the nwaitingnodes parameter be ignored at the root node?   */
   SCIP_Bool             uselprows;         /**< should subproblem be created out of the rows in the LP rows?      */
   SCIP_Bool             usegcg;            /**< should the subproblem be solved with GCG?                         */
   unsigned int          randseed;          /**< seed value for random number generator                            */
   SCIP_HASHTABLE*       hashtable;         /**< hashtable used to store the solution tuples already used          */
   SOLTUPLE*             lasttuple;         /**< last tuple of solutions created by crossover                      */
   int                   nrateredns;        /**< times similarityrate has already been reduced                     */
};

/** n-tuple of solutions and their hashkey */
struct SolTuple
{
   int*                  indices;            /**< sorted array of solution indices                                 */
   int                   size;               /**< size of the array (should be heurdata->nusedsols)                */
   unsigned int          key;                /**< hashkey of the tuple                                             */
   SOLTUPLE*             prev;               /**< previous solution tuple created                                  */
};




/*
 * Local methods
 */


/** get the points of which the current original solution is a convex combination */
static
SCIP_RETCODE getMembersOfDecomposition(
   SCIP*             scip,
   SCIP*             masterprob,
   SCIP_HEUR*        heur,
   int               npricingprobs,
   SCIP_SOL***       members,
   int*              nmembers
   )
{
   SCIP_VARDATA* vardata;
   SCIP_VARDATA* vardata2;
   int* memberind;
   int* blocknr;
   SCIP_Real* blockvalue;
   SCIP_Real increaseval;
   SCIP_VAR** mastervars;
   SCIP_Real* mastervals;
   int nmastervars;
   int i;
   int j;

   assert(scip != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &blockvalue, npricingprobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocknr, npricingprobs) );

   /* get variables of the master problem and their solution values */
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &mastervals, nmastervars) );
   SCIP_CALL( SCIPgetSolVals(masterprob, NULL, nmastervars, mastervars, mastervals) );

   SCIP_CALL( SCIPallocBufferArray(scip, &memberind, nmastervars) );

   /* initialize the block values for the pricing problems */
   for( i = 0; i < npricingprobs; i++ )
   {
      blockvalue[i] = 0.0;
      blocknr[i] = 0;
   }

   /* loop over all given master variables */
   for( i = 0; i < nmastervars; i++ )
   {
      memberind[i] = -1;

      /* first of all, handle the variables with integral values */
      while( SCIPisFeasGE(scip, mastervals[i], 1) )
      {
         vardata = SCIPvarGetData(mastervars[i]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_MASTER);
         assert(vardata->data.mastervardata.norigvars >= 0);
         assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
         assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);

         if( vardata->blocknr == -1 )
         {
            assert(vardata->data.mastervardata.norigvars == 2);
            assert(vardata->data.mastervardata.origvals[0] == 1.0);
            assert(vardata->data.mastervardata.origvals[1] == 0.0);

            (*nmembers)++;
            memberind[i] = *nmembers - 1;
            SCIP_CALL( SCIPreallocBufferArray(scip, members, *nmembers) );
            SCIP_CALL( SCIPcreateSol(scip, &(*members)[memberind[i]], heur) );

            /* increase the corresponding value */
            SCIP_CALL( SCIPincSolVal(scip, (*members)[memberind[i]], vardata->data.mastervardata.origvars[0], vardata->data.mastervardata.origvals[0]) );
            mastervals[i] = 0.0;
         }
         else
         {
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
            {
               if( SCIPisZero(scip, vardata->data.mastervardata.origvals[j]) )
                  continue;
               else if( memberind[i] == -1 )
               {
                  (*nmembers)++;
                  memberind[i] = *nmembers - 1;
                  SCIP_CALL( SCIPreallocBufferArray(scip, members, *nmembers) );
                  SCIP_CALL( SCIPcreateSol(scip, &(*members)[memberind[i]], heur) );
               }

               /* get the right original variable */
               vardata2 = SCIPvarGetData(vardata->data.mastervardata.origvars[j]);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);
               assert(vardata2->data.origvardata.pricingvar != NULL);
               vardata2 = SCIPvarGetData(vardata2->data.origvardata.pricingvar);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_PRICING);

               if( vardata2->data.pricingvardata.norigvars <= blocknr[vardata->blocknr] )
               {
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, (*members)[memberind[i]], vardata2->data.pricingvardata.origvars[vardata2->data.pricingvardata.norigvars-1], vardata->data.mastervardata.origvals[j]) );
                  mastervals[i] = 1.0;
               }
               else
               {
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, (*members)[memberind[i]], vardata2->data.pricingvardata.origvars[blocknr[vardata->blocknr]], vardata->data.mastervardata.origvals[j]) );
               }
            }
            mastervals[i] = mastervals[i] - 1.0;
            blocknr[vardata->blocknr]++;
         }
      }
   }

   /* loop over all given master variables */
   for( i = 0; i < nmastervars; i++ )
   {
      if( SCIPisFeasZero(scip, mastervals[i]) )
      {
         continue;
      }
      assert(SCIPisFeasGE(scip, mastervals[i], 0.0) && SCIPisFeasLT(scip, mastervals[i], 1.0));

      vardata = SCIPvarGetData(mastervars[i]);
      assert(vardata != NULL);
      assert(vardata->vartype == GCG_VARTYPE_MASTER);
      assert(vardata->data.mastervardata.norigvars >= 0);
      assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
      assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);

      while( SCIPisFeasPositive(scip, mastervals[i]) )
      {
         if( vardata->blocknr == -1 )
         {
            assert(vardata->data.mastervardata.norigvars == 2);
            assert(vardata->data.mastervardata.origvals[0] == 1.0);
            assert(vardata->data.mastervardata.origvals[1] == 0.0);

            if( memberind[i] == -1)
            {
               (*nmembers)++;
               SCIP_CALL( SCIPreallocBufferArray(scip, members, *nmembers) );
               SCIP_CALL( SCIPcreateSol(scip, &(*members)[*nmembers - 1], heur) );
               memberind[i] = *nmembers - 1;
            }

            /* increase the corresponding value */
            SCIP_CALL( SCIPincSolVal(scip, (*members)[memberind[i]], vardata->data.mastervardata.origvars[0], vardata->data.mastervardata.origvals[0]) );
            mastervals[i] = 0.0;
         }
         else
         {
            increaseval = MIN(mastervals[i], 1.0 - blockvalue[vardata->blocknr]);
            /* loop over all original variables contained in the current master variable */
            for( j = 0; j < vardata->data.mastervardata.norigvars; j++ )
            {
               if( SCIPisZero(scip, vardata->data.mastervardata.origvals[j]) )
                  continue;
               else if( memberind[i] == -1 )
               {
                  (*nmembers)++;
                  SCIP_CALL( SCIPreallocBufferArray(scip, members, *nmembers) );
                  SCIP_CALL( SCIPcreateSol(scip, &(*members)[*nmembers - 1], heur) );
                  memberind[i] = *nmembers - 1;
               }

               /* get the right original variable */
               vardata2 = SCIPvarGetData(vardata->data.mastervardata.origvars[j]);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_ORIGINAL);
               assert(vardata2->data.origvardata.pricingvar != NULL);
               vardata2 = SCIPvarGetData(vardata2->data.origvardata.pricingvar);
               assert(vardata2 != NULL);
               assert(vardata2->vartype == GCG_VARTYPE_PRICING);

               if( vardata2->data.pricingvardata.norigvars <= blocknr[vardata->blocknr] )
               {
                  increaseval = mastervals[i];

                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, (*members)[memberind[i]], vardata2->data.pricingvardata.origvars[vardata2->data.pricingvardata.norigvars-1], vardata->data.mastervardata.origvals[j]) );
               }
               else
               {
                  /* increase the corresponding value */
                  SCIP_CALL( SCIPincSolVal(scip, (*members)[memberind[i]], vardata2->data.pricingvardata.origvars[blocknr[vardata->blocknr]], vardata->data.mastervardata.origvals[j]) );
               }
            }

            mastervals[i] = mastervals[i] - increaseval;
            if( SCIPisFeasZero(scip, mastervals[i]) )
            {
               mastervals[i] = 0.0;
            }
            blockvalue[vardata->blocknr] += increaseval;

            /* if the value assigned to the block is equal to 1, this block is full and we take the next block */
            if( SCIPisFeasGE(scip, blockvalue[vardata->blocknr], 1.0) )
            {
               blockvalue[vardata->blocknr] = 0.0;
               blocknr[vardata->blocknr]++;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &mastervals);
   SCIPfreeBufferArray(scip, &memberind);
   SCIPfreeBufferArray(scip, &blocknr);
   SCIPfreeBufferArray(scip, &blockvalue);

   return SCIP_OKAY;
}

#if 0
/** prints members of decomposition to standard output */
static
SCIP_RETCODE printMembersOfDecomposition(
      SCIP*       scip,
      SCIP_SOL**  members,
      int         nmembers
      )
{
   int i;

//   printf("------------------------------------------------------------\n");
//   printf("Current LP solution:\n");
//   SCIP_CALL( SCIPprintSol(scip, GCGrelaxGetCurrentOrigSol(scip), NULL, FALSE) );
   printf("------------------------------------------------------------\n");
   for( i = 0; i < nmembers; i++)
   {
      printf("Member %i:\n", i);
      SCIP_CALL( SCIPprintSol(scip, members[i], NULL, FALSE) );
      printf("------------------------------------------------------------\n");
   }

   return SCIP_OKAY;
}
#endif

#if 0
/** sorts members of decomposition according to the objective */
static
SCIP_RETCODE sortMembersOfDecomposition(
   SCIP*          scip,
   SCIP_SOL**     members,
   int            nmembers
   )
{
   int i;
   int j;
   SCIP_SOL* tmp;

   /* simple insertion sort algorithm */
   for( i = 1; i < nmembers; i++ )
   {
      tmp = members[i];
      j = i-1;
      while( j >= 0
            && SCIPgetSolTransObj(scip, members[j]) > SCIPgetSolTransObj(scip, tmp) )
      {
         members[j+1] = members[j];
         j = j-1;
      }
      members[j+1] = tmp;
   }

   return SCIP_OKAY;
}
#else
/** sorts members of decomposition according to their feasibility
 ** TODO: This implementation is too expensive, the feasibilities may be calculated in advance */
static
SCIP_RETCODE sortMembersOfDecomposition(
   SCIP*          scip,
   SCIP_SOL**     members,
   int            nmembers
   )
{
   SCIP_ROW** rows;
   SCIP_Real* violations;
   int nrows;
   int feasibility;

   int i;
   int j;

   /* get LP rows of original problem */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   SCIP_CALL( SCIPallocBufferArray(scip, &violations, nmembers) );

   /* for each member, compute the total amount of row violations */
   for( i = 0; i < nmembers; i++ )
   {
      violations[i] = 0;
      for( j = 0; j < nrows; j++ )
      {
         feasibility = SCIPgetRowSolFeasibility(scip, rows[j], members[i]);
         if( SCIPisNegative(scip, feasibility) )
            violations[i] -= feasibility;
      }
   }

   SCIPsortRealPtr(violations, (void**) members, nmembers);

   SCIPfreeBufferArray(scip, &violations);

   return SCIP_OKAY;
}
#endif


/** local methods copied from crossover heuristic */

/** gets the hash key of a solution tuple */
static
SCIP_DECL_HASHGETKEY(hashGetKeySols)
{  /*lint --e{715}*/
   return elem;
}


/** returns TRUE iff both solution tuples are identical */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqSols)
{  /*lint --e{715}*/
   int i;
   int size;

   int* indices1;
   int* indices2;

   indices1 = ((SOLTUPLE*)key1)->indices;
   indices2 = ((SOLTUPLE*)key2)->indices;

   /* there should be two nonempty arrays of the same size */
   assert(indices1 != NULL);
   assert(indices2 != NULL);
   assert(((SOLTUPLE*)key1)->size == ((SOLTUPLE*)key2)->size);

   size  = ((SOLTUPLE*)key1)->size;

   /* compare arrays by components, return TRUE, iff equal */
   for( i = 0; i < size; i++ )
   {
      if( indices1[i] != indices2[i] )
         return FALSE;
   }

   return TRUE;
}


/** returns hashkey of a solution tuple */
static
SCIP_DECL_HASHKEYVAL(hashKeyValSols)
{  /*lint --e{715}*/
   return ((SOLTUPLE*)key)->key;
}


/** calculates a hash key for a given tuple of solution indices */
static
unsigned int calculateHashKey(
   int* indices,                            /**< indices of solutions                                            */
   int size                                 /**< number of solutions                                             */
   )
{
   int i;
   unsigned int hashkey;

   /* haskey should be (x1+1) * (x2+1) * ... * (xn+1) + x1 + x2 + ... + xn */
   hashkey = 1;
   for( i = 0; i < size; i++ )
      hashkey *= indices[i] + 1;
   for( i = 0; i < size; i++ )
      hashkey += indices[i];

   return hashkey;
}


/**< insertion sort for a small int array */
static void sortArray(
   int* a,                                  /**< array to be sorted                                              */
   int size                                 /**< size of array                                                   */
   )
{
   int i;
   int j;
   int tmp;

   /* simple insertion sort algorithm */
   for( i = 1; i < size; i++ )
   {
      tmp = a[i];
      j = i-1;
      while( j >= 0 && a[j] > tmp )
      {
         a[j+1] = a[j];
         j = j-1;
      }
      a[j+1] = tmp;
   }
}


/** creates a new tuple of solutions */
static
SCIP_RETCODE createSolTuple(
   SCIP*                 scip,              /**< original SCIP data structure                                    */
   SOLTUPLE**            elem,              /**< tuple of solutions which should be created                      */
   int*                  indices,           /**< indices of solutions                                            */
   int                   size,              /**< number of solutions                                             */
   SCIP_HEURDATA*        heurdata           /**< primal heuristic data                                           */
   )
{
   /* memory allociation */
   SCIP_CALL( SCIPallocBlockMemory(scip, elem) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*elem)->indices,size) );
   BMScopyMemoryArray((*elem)->indices, indices, size);

   /* data input */
   sortArray(indices,size);
   (*elem)->size = size;
   (*elem)->key = calculateHashKey((*elem)->indices, (*elem)->size);
   (*elem)->prev = heurdata->lasttuple;

   /* update heurdata */
   heurdata->lasttuple = *elem;
   return SCIP_OKAY;
}


/** randomly selects the solutions crossover will use from the pool of all solutions found so far */
static
SCIP_RETCODE selectSolsRandomized(
   SCIP*                 scip,              /**< original SCIP data structure                                    */
   int*                  selection,         /**< pool of solutions crossover uses                                */
   SCIP_HEURDATA*        heurdata,          /**< primal heuristic data                                           */
   int                   nmembers,          /**< number of members                                               */
   SCIP_Bool*            success            /**< pointer to store whether the proccess was successful            */
   )
{

   int i;
   int j;
   int lastsol;          /* the worst solution possible to choose */
   int nusedsols;        /* number of solutions which will be chosen */

   SOLTUPLE* elem;

   /* initialization */
   nusedsols = heurdata->nusedsols;
   lastsol = nmembers;
   assert(nusedsols < lastsol);

   i = 0;
   *success = FALSE;

   /* perform at maximum 100 restarts and stop as soon as a new set of solutions is found */
   while( !*success && i < 100 )
   {
      for( j = 0; j < nusedsols; j++ )
      {
         selection[j] = SCIPgetRandomInt(nusedsols-j-1, lastsol-1, &heurdata->randseed);
         lastsol = selection[j];
      }

      /* creates an object ready to be inserted into the hashtable */
      SCIP_CALL( createSolTuple(scip, &elem, selection, nusedsols, heurdata) );

      /* check whether the randomized set is already in the hashtable, if not, insert it */
      if( !SCIPhashtableExists(heurdata->hashtable, elem) )
      {
         SCIP_CALL( SCIPhashtableInsert(heurdata->hashtable, elem) );
         *success = TRUE;
      }
      i++;
   }

   return SCIP_OKAY;
}


/** check solutions for "similarity"; solutions are called similar if the rate of
 ** equal variables among all non-zero variables is >= heurdata->similarityrate */
static
SCIP_Bool isSimilar(
   SCIP*                scip,
   SCIP_HEURDATA*       heurdata,
   SCIP_SOL*            sol1,
   SCIP_SOL*            sol2
   )
{
   SCIP_VAR** vars;
   int nvars;

   int nequalvars;
   int nnonzerovars;
   int i;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   nequalvars = 0;
   nnonzerovars = 0;

   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real solval1;
      SCIP_Real solval2;

      solval1 = SCIPgetSolVal(scip, sol1, vars[i]);
      solval2 = SCIPgetSolVal(scip, sol2, vars[i]);

      if( !SCIPisZero(scip, solval1) || !SCIPisZero(scip, solval2) )
      {
         nnonzerovars++;
         if( SCIPisEQ(scip, solval1, solval2) )
            nequalvars++;
      }
   }

   if( nnonzerovars == 0 )
      return FALSE;

   if( SCIPisGE(scip, (SCIP_Real)nequalvars / (SCIP_Real)nnonzerovars, heurdata->similarityrate) )
      return TRUE;
   else
      return FALSE;
}

/** select solutions for crossover which are "similar" */
static
SCIP_RETCODE selectSimilarSols(
   SCIP*                scip,
   SCIP_HEURDATA*       heurdata,
   int*                 selection,
   SCIP_Bool*           used,
   SCIP_SOL**           members,
   int                  nmembers,
   SCIP_Bool*           success
   )
{
   int i;
   int j;
   int selsize;
   int nusedsols;

   nusedsols = heurdata->nusedsols;

   *success = FALSE;

   do
   {
      selsize = 0;
      /* loop over all members and break if we find an acceptable solution set for crossover;
       * members which were tried to be selected before are ignored */
      for( i = 0; i < nmembers - 1 && !(*success); i++ )
      {
         SOLTUPLE* elem;

         if( used[i] )
            continue;

//         SCIPdebugMessage("Trying to find similar solutions:\n");
         selection[0] = i;
         used[i] = TRUE;
         selsize = 1;
//         SCIPdebugMessage("  -> member %i.\n", i);

         /* loop over all remaining members and try to add the similar ones to the selection */
         for( j = i + 1; j < nmembers && selsize < nusedsols; j++ )
         {
            if( isSimilar(scip, heurdata, members[i], members[j]) && !used[j] )
            {
               selection[selsize] = j;
               used[j] = TRUE;
               selsize++;
//               SCIPdebugMessage("  -> member %i.\n", j);
            }
         }

         /* if the selection is large enough for crossover, check if it has already been used before */
         if( selsize == nusedsols )
         {
            SCIP_CALL( createSolTuple(scip, &elem, selection, nusedsols, heurdata) );
            if( !SCIPhashtableExists(heurdata->hashtable, elem) )
            {
               SCIP_CALL( SCIPhashtableInsert(heurdata->hashtable, elem) );
               *success = TRUE;
            }
            else
               continue;
         }
      }

      /* if no solution set for crossover could be found, try to reduce similarityrate and start again */
      if( !(*success) && heurdata->nrateredns < heurdata->maxrateredns )
      {
         for( i = 0; i < nmembers; i++ )
            used[i] = FALSE;
         heurdata->similarityrate *= 0.75;
         heurdata->nrateredns++;
//         SCIPdebugMessage("no solution set for crossover found - reducing similarityrate to %.2f\n", heurdata->similarityrate);
      }
   }
   while( !(*success) && heurdata->nrateredns < heurdata->maxrateredns );

   assert(selsize <= nusedsols);

   if( !(*success) )
   {
      if( heurdata->randomization )
      {
//         SCIPdebugMessage("No similar solutions found - selecting randomly\n");
         SCIP_CALL( selectSolsRandomized(scip, selection, heurdata, nmembers, success) );
      }
      else
      {
//         SCIPdebugMessage("No similar solutions found.\n");
      }
   }

   return SCIP_OKAY;
}


/** creates the all variables of the subproblem */
static SCIP_RETCODE fixVariables(
   SCIP*                 scip,               /**< original SCIP data structure                                  */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                        */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                               */
   int*                  selection,          /**< pool of solutions crossover will use                          */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data                                         */
   SCIP_SOL**            members,            /**< members of decomposition                                      */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_VAR** vars;                          /* original scip variables                */
   SCIP_Real fixingrate;                     /* percentage of variables that are fixed */

   int nvars;
   int nbinvars;
   int nintvars;

   int i;
   int j;
   int fixingcounter;
   int zerocounter;

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   fixingcounter = 0;
   zerocounter = 0;

   /* create the binary and general integer variables of the subproblem */
   for( i = 0; i < nbinvars + nintvars; i++ )
   {
      SCIP_Real solval;
      SCIP_Real lpsolval;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Bool fixable;

      fixable = TRUE;
      solval = SCIPgetSolVal(scip, members[selection[0]], vars[i]);

      /* check, whether variable's value is identical for each selected solution */
      for( j = 1; j < heurdata->nusedsols; j++ )
      {
         SCIP_Real varsolval;
         varsolval = SCIPgetSolVal(scip, members[selection[j]], vars[i]);
         if( REALABS(solval - varsolval) > 0.5 )
         {
            fixable = FALSE;
            break;
         }
      }

//      /* if fixing the variable would exceed maxfixingrate, don't fix it */
//      if( (SCIP_Real)(fixingcounter + 1) / (SCIP_Real)(MAX(nbinvars + nintvars, 1)) > heurdata->maxfixingrate )
//         fixable = FALSE;

      /* if variable is 0 for all members in selection, but != 0 for LP optimum, don't fix it */
      lpsolval = SCIPgetSolVal(scip, GCGrelaxGetCurrentOrigSol(scip), vars[i]);
      if( fixable && SCIPisZero(scip, solval) && !SCIPisZero(scip, lpsolval) )
         fixable = FALSE;

      /* fix other zero variables only with a certain probability */
      if( fixable && SCIPisZero(scip, solval) )
      {
         if( SCIPgetRandomReal(0, 1, &heurdata->randseed) >= heurdata->zerofixprob )
            fixable = FALSE;
      }

      /* get variable's global bounds */
      lb = SCIPvarGetLbGlobal(vars[i]);
      ub = SCIPvarGetUbGlobal(vars[i]);

      /* original solval can be outside transformed global bounds */
      if( fixable && (lb > solval || solval > ub) )
         fixable = FALSE;

      /* if solutions' values are equal, variable is fixed in the subproblem, otherwise it is just copied */
      if( fixable )
      {
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], solval) );
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], solval) );
         fixingcounter++;

         if ( SCIPisZero(scip, solval) )
            zerocounter++;
      }
   }

   fixingrate = (SCIP_Real)fixingcounter / (SCIP_Real)(MAX(nbinvars + nintvars, 1));

   SCIPdebugMessage("subSCIP: %i out of %i (%.2f percent) variables have been fixed.\n", fixingcounter, nbinvars + nintvars, fixingrate * 100.0);
   SCIPdebugMessage("subSCIP: %i out of %i (%.2f percent) fixed variables are zero.\n", zerocounter, fixingcounter,
         (SCIP_Real)zerocounter / (SCIP_Real)fixingcounter * 100.0);

   /* if all variables were fixed or amount of fixed variables is insufficient, skip residual part of
    * subproblem creation ans abort immediately */
   if( fixingcounter == nbinvars + nintvars || fixingrate < heurdata->minfixingrate )
   {
      *success = FALSE;
      SCIPdebugMessage("Fixing of variables was not successful - fixing rate %f percent.\n", fixingrate * 100.0);
   }

   return SCIP_OKAY;
}

/** creates the rows of the subproblem */
static
SCIP_RETCODE createRows(
   SCIP*                 scip,               /**< original SCIP data structure                                  */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                        */
   SCIP_VAR**            subvars             /**< the variables of the subproblem                               */
   )
{
   SCIP_ROW** rows;                          /* original scip rows                       */
   SCIP_CONS* cons;                          /* new constraint                           */
   SCIP_VAR** consvars;                      /* new constraint's variables               */
   SCIP_COL** cols;                          /* original row's columns                   */

   SCIP_Real constant;                       /* constant added to the row                */
   SCIP_Real lhs;                            /* left hand side of the row                */
   SCIP_Real rhs;                            /* left right side of the row               */
   SCIP_Real* vals;                          /* variables' coefficient values of the row */

   int nrows;
   int nnonz;
   int i;
   int j;

   /* get the rows and their number */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* copy all rows to linear constraints */
   for( i = 0; i < nrows; i++ )
   {
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

      assert(lhs <= rhs);

      /* allocate memory array to be filled with the corresponding subproblem variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nnonz) );
      for( j = 0; j < nnonz; j++ )
         consvars[j] = subvars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

      /* create a new linear constraint and add it to the subproblem */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );

      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** creates a subproblem for subscip by fixing a number of variables */
static
SCIP_RETCODE setupSubproblem(
   SCIP*                 scip,               /**< original SCIP data structure                                  */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                        */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                               */
   int*                  selection,          /**< pool of solutions crossover will use                          */
   SCIP_Bool*            used,               /**< vertices that have already been used for crossover            */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data                                         */
   SCIP_SOL**            members,            /**< members of decomposition                                      */
   int                   nmembers,           /**< number of members                                             */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{

   int nusedsols;                           /* number of solutions to use in crossover     */

//   int i;

   /* get solutions' data */
   nusedsols = heurdata->nusedsols;

   assert(nusedsols > 1);
   assert(nmembers >= nusedsols);

#if 0
   /* use nusedsols best solutions if randomization is deactivated or there are only nusedsols solutions at hand
    * or a good new solution was found since last call */
   if( !heurdata->randomization || nmembers == nusedsols || heurdata->prevlastsol != members[nusedsols-1] )
   {
      SOLTUPLE* elem;

      for( i = 0; i < nusedsols; i++ )
         selection[i] = i;
      SCIP_CALL( createSolTuple(scip, &elem, selection, nusedsols, heurdata) );

      /* check, whether solution tuple has already been tried */
      if( SCIPhashtableExists(heurdata->hashtable, elem) )
      {
         *success = FALSE;

         /* if solution tuple has already been tried, randomization is allowed and enough solutions are at hand, try
          * to randomize another tuple. E.g., this can happen if the last crossover solution was among the best ones */
         if( heurdata->randomization && nmembers > nusedsols )
         {
            SCIP_CALL( selectSolsRandomized(scip, selection, heurdata, nmembers, success) );
         }
      }
      else
      {
         *success = TRUE;
         SCIP_CALL( SCIPhashtableInsert(heurdata->hashtable, elem) );
      }
   }
   /* otherwise randomize the set of solutions */
   else
   {
      SCIP_CALL( selectSolsRandomized(scip, selection, heurdata, nmembers, success) );
   }
#else
   SCIP_CALL( selectSimilarSols(scip, heurdata, selection, used, members, nmembers, success) );
#endif

   /* no acceptable solution tuple could be created */
   if( !(*success) )
      return SCIP_OKAY;

   /* create the variables of the subproblem */
   SCIP_CALL( fixVariables(scip, subscip, subvars, selection, heurdata, members, success) );

   /* if the enough variables could be fixed, create rows of the subproblem */
   if( *success && heurdata->uselprows )
   {
      SCIP_CALL( createRows(scip, subscip, subvars) );
   }

   return SCIP_OKAY;
}


/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< crossover heuristic structure                       */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   int*                  solindex,           /**< index of the solution                               */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );
   *solindex = SCIPsolGetIndex(newsol);

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, success) );

   if( *success )
   {
      SCIPdebugMessage("GCG extreme points crossover: new solution added.\n");
   }

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** updates heurdata after a run of crossover */
static
void updateFailureStatistic(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data                               */
   )
{
   /* increase number of failures, calculate next node at which crossover should be called and update actual solutions */
   heurdata->nfailures++;
   heurdata->nextnodenumber = (heurdata->nfailures <= 25
      ? SCIPgetNNodes(scip) + 100*(2LL << heurdata->nfailures)
      : SCIP_LONGINT_MAX);
}

/** applies crossover on members of decomposition */
static
SCIP_RETCODE applyCrossover(
      SCIP*                scip,
      SCIP_HEUR*           heur,
      SCIP_RESULT*         result,
      SCIP_Real            memorylimit,
      SCIP_Real            timelimit,
      SCIP_Longint         nstallnodes,
      SCIP_SOL**           members,
      int                  nmembers
      )
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** vars;                          /* original problem's variables                        */

   SCIP_Real cutoff;                         /* objective cutoff for the subproblem                 */
   SCIP_Real upperbound;
   SCIP_Bool success;

   SCIP_Bool* used;
   int nvars;                                /* number of original problem's variables              */
   int nusedsols;
   int maxruns;
//   int maxsize;

   int i;
   int runs;

#ifdef NDEBUG
   SCIP_RETCODE retstat;
#endif

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   nusedsols = heurdata->nusedsols;
   maxruns = heurdata->maxruns;

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPallocBufferArray(scip, &used, nmembers) );
   for( i = 0; i < nmembers; i++ )
      used[i] = FALSE;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(nvars > 0);

   success = FALSE;

   for( runs = 0; runs < maxruns && !success; runs++ )
   {

      SCIP* subscip;                            /* the subproblem created by crossover                 */
      SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to subSCIP variables      */
      SCIP_VAR** subvars;                       /* subproblem's variables                              */
      int* selection;                           /* pool of solutions crossover uses                    */

      SCIPdebugMessage("Crossover -- run %i (minfixingrate = %.2f, similarityrate = %.2f)\n",
            runs + 1, heurdata->minfixingrate, heurdata->similarityrate);

      /* initialize the subproblem */
      SCIP_CALL( SCIPcreate(&subscip) );

      SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &selection, nusedsols) );

      if( heurdata->usegcg )
      {
         if( heurdata->uselprows )
         {
            char probname[SCIP_MAXSTRLEN];

            /* copy all plugins */
            SCIP_CALL( SCIPincludeGcgPlugins(subscip) );

            /* get name of the original problem and add the string "_extremeptsub" */
            (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_extremeptsub", SCIPgetProbName(scip));

            /* create the subproblem */
            SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

            /* copy all variables */
            SCIP_CALL( SCIPcopyVars(scip, subscip, varmapfw, NULL, TRUE) );
         }
         else
         {
            SCIP_Bool valid;
            valid = FALSE;

            SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "extremept", TRUE, FALSE, &valid) );
            SCIPdebugMessage("Copying the SCIP instance was %s complete.\n", valid ? "" : "not ");
         }

         // TODO: copy the matrix structure here
         *result = SCIP_DIDNOTFIND;
         return SCIP_OKAY;
      }
      else
      {
         char probname[SCIP_MAXSTRLEN];

         /* copy all plugins */
         SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

         /* get name of the original problem and add the string "_extremeptsub" */
         (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_extremeptsub", SCIPgetProbName(scip));

         /* create the subproblem */
         SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

         /* copy all variables */
         SCIP_CALL( SCIPcopyVars(scip, subscip, varmapfw, NULL, TRUE) );

         /* if the lp rows are not used, also copy the constraints */
         if( !heurdata->uselprows )
         {
            SCIP_Bool valid;
            valid = FALSE;

            SCIP_CALL( SCIPcopyConss(scip, subscip, varmapfw, NULL, TRUE, FALSE, &valid) );
            SCIPdebugMessage("Copying the SCIP constraints was %s complete.\n", valid ? "" : "not ");
         }
      }

      for( i = 0; i < nvars; i++ )
         subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

      /* free hash map */
      SCIPhashmapFree(&varmapfw);

      /* create a new problem, which fixes variables with same value in bestsol and LP relaxation */
      SCIP_CALL( setupSubproblem(scip, subscip, subvars, selection, used, heurdata, members, nmembers, &success) );
      SCIPdebugMessage("Extreme Points Crossover subproblem: %d vars, %d cons, success=%u\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip), success);

      /* if creation of subscip was aborted (e.g. due to number of fixings), free subscip and abort */
      if( !success )
      {
//         *result = SCIP_DIDNOTRUN;

         /* if creation was aborted due to number of fixings, free the already created subproblem */
         if( SCIPgetStage(subscip) != SCIP_STAGE_INIT )
         {
            int nbinvars;
            int nintvars;
            SCIP_CALL( SCIPgetVarsData(scip, NULL, NULL, &nbinvars, &nintvars, NULL, NULL) );
            for( i = 0; i < nbinvars + nintvars; i++ )
            {
               SCIP_CALL( SCIPreleaseVar(subscip, &subvars[i]) );
            }
            SCIP_CALL( SCIPfreeTransform(subscip) );
         }

         /* free memory */
         SCIPfreeBufferArray(scip, &selection);
         SCIPfreeBufferArray(scip, &subvars);
         SCIP_CALL( SCIPfree(&subscip) );

         /* this run will be counted as a failure since no new solution tuple could be generated or the neighborhood of the
          * solution was not fruitful in the sense that it was too big
          */
         updateFailureStatistic(scip, heurdata);

         return SCIP_OKAY;
      }

      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

      /* disable output to console */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

      /* set limits for the subproblem */
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nstallnodes) );
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

      /* forbid recursive call of heuristics and separators solving subMIPs */
      SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

      /* disable cutting plane separation */
      SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

      /* disable expensive presolving */
      SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

      /* use best estimate node selection */
      if( SCIPfindNodesel(scip, "estimate") != NULL )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
      }

      /* use inference branching */
      if( SCIPfindBranchrule(scip, "inference") != NULL )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
      }

      /* disable conflict analysis */
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useprop", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useinflp", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useboundlp", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usesb", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usepseudo", FALSE) );

      /* if there is already a solution, add an objective cutoff */
      if( SCIPgetNSols(scip) > 0 )
      {
         cutoff = SCIPinfinity(scip);
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
         cutoff = MIN(upperbound, cutoff );
         SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );
      }

      /* solve the subproblem */
      SCIPdebugMessage("subSCIP: Solving... (node limit = %lld, time limit = %.2g)\n", nstallnodes, timelimit);

      /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
       * Hence in optimized mode, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
       */
#ifdef NDEBUG
      retstat = SCIPsolve(subscip);
      if( retstat != SCIP_OKAY )
      {
         SCIPwarningMessage("Error while solving subMIP in GCG extreme points crossover heuristic; subSCIP terminated with code <%d>\n",
               retstat);
      }
#else
      SCIP_CALL( SCIPsolve(subscip) );
#endif

      heurdata->usednodes += SCIPgetNNodes(subscip);

      /* check, whether a solution was found */
      success = FALSE;
      if( SCIPgetNSols(subscip) > 0 )
      {
         SCIP_SOL** subsols;
         int nsubsols;
         int solindex;                             /* index of the solution created by crossover          */

         SCIPdebugMessage("GCG extreme points crossover heuristic found %i feasible solution(s).\n", SCIPgetNSols(subscip));

         /* check, whether a solution was found;
          * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
          */
         nsubsols = SCIPgetNSols(subscip);
         subsols = SCIPgetSols(subscip);
         solindex = -1;
         for( i = 0; i < nsubsols && !success; ++i )
         {
            SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &solindex, &success) );
         }

         if( success )
            *result = SCIP_FOUNDSOL;
         else
            updateFailureStatistic(scip, heurdata);
      }
      else
      {
         /* if no new solution was found, run was a failure */
         updateFailureStatistic(scip, heurdata);
         SCIPdebugMessage("GCG extreme points crossover: no subMIP solution found - ");
         switch ( SCIPgetStatus(subscip) ) {
         case SCIP_STATUS_INFEASIBLE:
            SCIPdebugPrintf("subMIP infeasible.\n");
            heurdata->minfixingrate *= 0.8;
            break;
         case SCIP_STATUS_NODELIMIT:
         case SCIP_STATUS_STALLNODELIMIT:
            SCIPdebugPrintf("node limit reached.\n");
            heurdata->minfixingrate *= 1.2;
            break;
         case SCIP_STATUS_TIMELIMIT:
            SCIPdebugPrintf("time limit reached.\n");
            heurdata->minfixingrate *= 1.2;
            break;
         case SCIP_STATUS_USERINTERRUPT:
            SCIPdebugPrintf("solving process interrupted by user.\n");
            break;
         default:
            SCIPdebugPrintf("SCIP status %d.\n", SCIPgetStatus(subscip));
            break;
         }
      }

      /* free subproblem */
      SCIPfreeBufferArray(scip, &selection);
      SCIP_CALL( SCIPfreeTransform(subscip) );
      for( i = 0; i < nvars; i++ )
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &subvars[i]) );
      }
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
   }

   SCIPfreeBufferArray(scip, &used);

   return SCIP_OKAY;
}



/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#define heurCopyExtremepoints NULL

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeExtremepoints)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitExtremepoints)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->randseed = 0;
   heurdata->lasttuple = NULL;
   heurdata->nfailures = 0;
   heurdata->nextnodenumber = 0;
   heurdata->nrateredns = 0;

   /* initialize hash table */
   SCIP_CALL( SCIPhashtableCreate(&heurdata->hashtable, SCIPblkmem(scip), HASHSIZE_SOLS,
         hashGetKeySols, hashKeyEqSols, hashKeyValSols, NULL) );
   assert(heurdata->hashtable != NULL );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitExtremepoints)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SOLTUPLE* soltuple;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   soltuple = heurdata->lasttuple;

   /* free all soltuples iteratively */
   while( soltuple != NULL )
   {
      SOLTUPLE* tmp;
      tmp = soltuple->prev;
      SCIPfreeBlockMemoryArray(scip,&soltuple->indices,soltuple->size);
      SCIPfreeBlockMemory(scip,&soltuple);
      soltuple = tmp;
   }

   /* free hash table */
   assert(heurdata->hashtable != NULL );
   SCIPhashtableFree(&heurdata->hashtable);

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolExtremepoints NULL


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolExtremepoints NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecExtremepoints)
{  /*lint --e{715}*/

   SCIP* masterprob;
   SCIP_HEURDATA* heurdata;

   SCIP_Real memorylimit;                    /* memory limit for the subproblem                     */
   SCIP_Real timelimit;                      /* time limit for the subproblem                       */
   SCIP_Longint nstallnodes;                 /* node limit for the subproblem                       */

   int nmembers;
   SCIP_SOL** members;                       /* array to store relevant extreme points              */

   int i;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /* get master problem */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   assert(SCIPhasCurrentNodeLP(masterprob));

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result= SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(masterprob) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   *result = SCIP_DELAYED;

   /* if heuristic should be delayed, wait until certain number of nodes is reached */
   if( SCIPgetNNodes(scip) < heurdata->nextnodenumber )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward Crossover if it succeeded often */
   nstallnodes = (SCIP_Longint)
                              (nstallnodes * (1.0 + 2.0*(SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur)+1.0)));

   /* count the setup costs for the sub-MIP as 100 nodes */
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMessage("skipping GCG extreme points crossover: nstallnodes=%"SCIP_LONGINT_FORMAT", minnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if( timelimit < 10.0 || memorylimit <= 0.0 )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPdebugMessage("Executing GCG extreme points crossover heuristic ...\n");

   /* get members of decomposition, i. e. the points of which the relaxation solution is a convex combination */

   SCIP_CALL( SCIPallocBufferArray(scip, &members, 1) );
   nmembers = 0;
   SCIP_CALL( getMembersOfDecomposition(scip, masterprob, heur, GCGrelaxGetNPricingprobs(scip), &members, &nmembers) );
   SCIPdebugMessage("Current original LP solution consists of %i members of decomposition.\n", nmembers);

   *result = SCIP_DIDNOTFIND;

   /* only call heuristic if there are enough members */
   if( nmembers < heurdata->nusedsols  ) {
      /* free memory */
      for ( i = 0; i < nmembers; i++ )
      {
         SCIP_CALL( SCIPfreeSol(scip, &members[i]) );
      }
      SCIPfreeBufferArray(scip, &members);
      return SCIP_OKAY;
   }

   SCIP_CALL( sortMembersOfDecomposition(scip, members, nmembers) );

//#ifdef SCIP_DEBUG
//   SCIP_CALL( printMembersOfDecomposition(scip, members, nmembers) );
//#endif

   SCIP_CALL( applyCrossover(scip, heur, result, memorylimit, timelimit, nstallnodes, members, nmembers) );

   /* free memory */
   for ( i = 0; i < nmembers; i++ )
   {
      SCIP_CALL( SCIPfreeSol(scip, &members[i]) );
   }
   SCIPfreeBufferArray(scip, &members);

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the extreme points crossover primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurExtremepoints(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create extreme points crossover primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyExtremepoints, heurFreeExtremepoints, heurInitExtremepoints, heurExitExtremepoints,
         heurInitsolExtremepoints, heurExitsolExtremepoints, heurExecExtremepoints,
         heurdata) );

   /* add extreme points crossover primal heuristic parameters */

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxruns",
         "maximum times the heuristic runs on an LP solution",
         &heurdata->maxruns, FALSE, DEFAULT_MAXRUNS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/nusedsols",
         "number of solutions to be taken into account",
         &heurdata->nusedsols, FALSE, DEFAULT_NUSEDSOLS, 2, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/similarityrate",
         "rate by which to decide whether two solutions are \"similar\"",
         &heurdata->similarityrate, FALSE, DEFAULT_SIMILARITYRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxrateredns",
         "maximum times that similarityrate is allowed to be reduced",
         &heurdata->maxrateredns, TRUE, DEFAULT_MAXRATEREDNS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minfixingrate",
         "minimum percentage of integer variables that have to be fixed",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/maxfixingrate",
         "maximum percentage of integer variables that can be fixed",
         &heurdata->maxfixingrate, FALSE, DEFAULT_MAXFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/zerofixprob",
         "probability that a zero variable is fixed",
         &heurdata->zerofixprob, FALSE, DEFAULT_ZEROFIXPROB, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which crossover should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/randomization",
         "should the choice which sols to take be randomized?",
         &heurdata->randomization, TRUE, DEFAULT_RANDOMIZATION, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/dontwaitatroot",
         "should the nwaitingnodes parameter be ignored at the root node?",
         &heurdata->dontwaitatroot, TRUE, DEFAULT_DONTWAITATROOT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/usegcg",
         "should the subproblem be solved with GCG?",
         &heurdata->usegcg, FALSE, DEFAULT_USEGCG, NULL, NULL) );

   return SCIP_OKAY;
}
