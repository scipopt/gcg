/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_master.c
 * @brief  master separator
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <stdio.h>

#include "scip/scip.h"
#include "scip/lp.h"
#include "sepa_master.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"


#define SEPA_NAME         "master"
#define SEPA_DESC         "separator for separating cuts in the original problem, called in the master"
#define SEPA_PRIORITY     0

#define SEPA_FREQ         1
#define SEPA_MAXBOUNDDIST 1.0
#define SEPA_USESSUBSCIP  FALSE
#define SEPA_DELAY        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define STARTMAXCUTS 50       /**< maximal cuts used at the beginning */
#define MAXCUTSINC   20       /**< increase of allowed number of cuts */


/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_ROW** mastercuts;     /**< cuts in the master problem */
   SCIP_ROW** origcuts;       /**< cuts in the original problem */
   int        norigcuts;      /**< number of cuts in the original problem */
   int        nmastercuts;    /**< number of cuts in the master problem */
   int        maxcuts;        /**< maximal number of allowed cuts */
};


/*
 * Local methods
 */


/** allocates enough memory to hold more cuts */
static
SCIP_RETCODE ensureSizeCuts(
   SCIP*          scip,     /**< SCIP data structure */  
   SCIP_SEPADATA* sepadata, /**< separator data data structure */
   int            size      /**< new size of cut arrays */  
   )
{
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(sepadata->mastercuts != NULL);
   assert(sepadata->origcuts != NULL);
   assert(sepadata->norigcuts <= sepadata->maxcuts);
   assert(sepadata->norigcuts >= 0);
   assert(sepadata->nmastercuts <= sepadata->maxcuts);
   assert(sepadata->nmastercuts >= 0);

   if( sepadata->maxcuts < size )
   {
      while ( sepadata->maxcuts < size )
      {
         sepadata->maxcuts += MAXCUTSINC;
      }
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(sepadata->mastercuts), sepadata->maxcuts) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(sepadata->origcuts), sepadata->maxcuts) );
   }
   assert(sepadata->maxcuts >= size);

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

#define sepaCopyMaster NULL
#define sepaInitMaster NULL
#define sepaInitsolMaster NULL
#define sepaExecsolMaster NULL

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeMaster)
{
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);

   SCIPfreeMemoryArray(scip, &(sepadata->origcuts));
   SCIPfreeMemoryArray(scip, &(sepadata->mastercuts));

   SCIPfreeMemory(scip, &sepadata);

   return SCIP_OKAY;
}

/** deinitialization method of separator (called before transformed problem is freed) */
static
SCIP_DECL_SEPAEXIT(sepaExitMaster)
{
   SCIP* origscip;
   SCIP_SEPADATA* sepadata;
   int i;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->nmastercuts == sepadata->norigcuts);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   for( i = 0; i < sepadata->norigcuts; i++ )
   {
      SCIP_CALL( SCIPreleaseRow(origscip, &(sepadata->origcuts[i])) );
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolMaster)
{
   SCIP_SEPADATA* sepadata;
   int i;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->nmastercuts == sepadata->norigcuts);

   assert(GCGpricerGetOrigprob(scip) != NULL);

   for( i = 0; i < sepadata->nmastercuts; i++ )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(sepadata->mastercuts[i])) );
   }

   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpMaster)
{
   SCIP*   origscip;
   SCIP_Bool delayed;
   SCIP_Bool cutoff;
   SCIP_ROW** cuts;
   SCIP_ROW* mastercut;
   SCIP_ROW* origcut;
   SCIP_COL** cols;
   SCIP_VAR** rowvars;
   int ncols;
   SCIP_Real* vals;
   SCIP_SEPADATA* sepadata;

   SCIP_VAR** mastervars;
   SCIP_Real* mastervals;
   int nmastervars;

   char name[SCIP_MAXSTRLEN];

   int ncuts;
   int i;
   int j;
   SCIP_Bool feasible;

   assert(scip != NULL);
   assert(result != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPdebugMessage("sepaExeclpMaster\n");

   *result = SCIP_DIDNOTFIND;

   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("master LP not solved to optimality, do no separation!\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( GCGrelaxUpdateCurrentSol(origscip, &feasible) );

   SCIP_CALL( SCIPseparateSol(origscip, GCGrelaxGetCurrentOrigSol(origscip),
         FALSE, FALSE, &delayed, &cutoff) );

   SCIPdebugMessage("SCIPseparateSol() found %d cuts!\n", SCIPgetNCuts(origscip));

   cuts = SCIPgetCuts(origscip);
   ncuts = SCIPgetNCuts(origscip);

   /* save cuts in the origcuts array in the separator data */
   assert(sepadata->norigcuts == sepadata->nmastercuts);
   SCIP_CALL( ensureSizeCuts(scip, sepadata, sepadata->norigcuts + ncuts) );
   for( i = 0; i < ncuts; i++ )
   {
      sepadata->origcuts[sepadata->norigcuts] = cuts[i];
      SCIP_CALL( SCIPcaptureRow(origscip, sepadata->origcuts[sepadata->norigcuts]) );
      sepadata->norigcuts++;
   }

   SCIP_CALL( SCIPclearCuts(origscip) );

   mastervars = SCIPgetVars(scip);
   nmastervars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &mastervals, nmastervars) );

   for( i = 0; i < ncuts; i++ )
   {
      origcut = sepadata->origcuts[sepadata->norigcuts-ncuts+i]; /*lint !e679*/ 

      /* get columns and vals of the cut */
      ncols = SCIProwGetNNonz(origcut);
      cols = SCIProwGetCols(origcut);
      vals = SCIProwGetVals(origcut);

      /* get the variables corresponding to the columns in the cut */
      SCIP_CALL( SCIPallocBufferArray(scip, &rowvars, ncols) );
      for( j = 0; j < ncols; j++ )
      {
         rowvars[j] = SCIPcolGetVar(cols[j]);
      }

      /* create new cut in the master problem */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "mc_%s", SCIProwGetName(origcut));
      SCIP_CALL( SCIPcreateEmptyRow(scip, &mastercut, name,
            ( SCIPisInfinity(scip, -SCIProwGetLhs(origcut)) ?
               SCIProwGetLhs(origcut) : SCIProwGetLhs(origcut) - SCIProwGetConstant(origcut)),
            ( SCIPisInfinity(scip, SCIProwGetRhs(origcut)) ?
               SCIProwGetRhs(origcut) : SCIProwGetRhs(origcut) - SCIProwGetConstant(origcut)),
            SCIProwIsLocal(origcut), TRUE, FALSE) );

      /* transform the original variables to master variables and add them to the cut */
      GCGrelaxTransformOrigvalsToMastervals(GCGpricerGetOrigprob(scip), rowvars, vals, ncols, mastervars, mastervals, nmastervars);
      SCIP_CALL( SCIPaddVarsToRow(scip, mastercut, nmastervars, mastervars, mastervals) );

      /* add the cut to the master problem */
      SCIP_CALL( SCIPaddCut(scip, NULL, mastercut, FALSE) );
      sepadata->mastercuts[sepadata->nmastercuts] = mastercut;
      sepadata->nmastercuts++;

#ifdef SCIP_DEBUG
      //SCIPdebugMessage("Cut %d:\n", i);
      //SCIP_CALL( SCIPprintRow(scip, mastercut, NULL) );
      //SCIPdebugMessage("\n\n");
#endif

      SCIPfreeBufferArray(scip, &rowvars);
   }

   if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   SCIPdebugMessage("%d cuts are in the original sepastore!\n", SCIPgetNCuts(origscip));
   SCIPdebugMessage("%d cuts are in the master sepastore!\n", SCIPgetNCuts(scip));

   SCIPfreeBufferArray(scip, &mastervals);

   assert(sepadata->norigcuts == sepadata->nmastercuts );

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the master separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaMaster(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create master separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(sepadata->origcuts), STARTMAXCUTS) ); /*lint !e506*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sepadata->mastercuts), STARTMAXCUTS) ); /*lint !e506*/
   sepadata->maxcuts = STARTMAXCUTS;
   sepadata->norigcuts = 0;
   sepadata->nmastercuts = 0;

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaCopyMaster, sepaFreeMaster, sepaInitMaster, sepaExitMaster,
         sepaInitsolMaster, sepaExitsolMaster,
         sepaExeclpMaster, sepaExecsolMaster,
         sepadata) );

   return SCIP_OKAY;
}


/* returns the array of original cuts saved in the separator data */
SCIP_ROW** GCGsepaGetOrigcuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);

   sepa = SCIPfindSepa(scip, SEPA_NAME);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   return sepadata->origcuts;
}

/* returns the number of original cuts saved in the separator data */
int GCGsepaGetNOrigcuts(
   SCIP* scip                 /**< SCIP data structure */
   )
{
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);

   sepa = SCIPfindSepa(scip, SEPA_NAME);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   return sepadata->norigcuts;
}

/* returns the array of master cuts saved in the separator data */
SCIP_ROW** GCGsepaGetMastercuts(
   SCIP* scip                 /**< SCIP data structure */
   )
{
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);

   sepa = SCIPfindSepa(scip, SEPA_NAME);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   return sepadata->mastercuts;
}

/* returns the number of master cuts saved in the separator data */
int GCGsepaGetNMastercuts(
   SCIP* scip                 /**< SCIP data structure */
   )
{
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);

   sepa = SCIPfindSepa(scip, SEPA_NAME);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   return sepadata->nmastercuts;
}
