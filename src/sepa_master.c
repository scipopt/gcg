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
/**@file   sepa_master.c
 * @ingroup SEPARATORS
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


#define SEPA_NAME              "master"
#define SEPA_DESC              "separator for separating cuts in the original problem, called in the master"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                     1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP            FALSE
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */


#define STARTMAXCUTS 50
#define MAXCUTSINC 20


/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_ROW** mastercuts;
   SCIP_ROW** origcuts;
   int norigcuts;
   int nmastercuts;
   int maxcuts;  
};




/*
 * Local methods
 */

static
SCIP_RETCODE ensureSizeCuts(
   SCIP*                 scip,
   SCIP_SEPADATA*        sepadata,
   int                   size
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

#if 0
static
SCIP_RETCODE checkCutConsistency(
   SCIP*                 scip
   )
{
   SCIP_SEPA* sepa;
   SCIP_SEPADATA* sepadata;
   int i;
   int j;
   int k;
   int l;
   int v;
   SCIP_COL** origcols;
   SCIP_COL** mastercols;
   SCIP_Real* origvals;
   SCIP_Real* mastervals;
   int norigcols;
   int nmastercols;
   
   SCIP_Real sum;
   SCIP_Real masterval;

   SCIP_VAR** mastervars;
   int nmastervars;
   SCIP_VARDATA* vardata;

   sepa = SCIPfindSepa(scip, SEPA_NAME);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   mastervars = SCIPgetVars(scip);
   nmastervars = SCIPgetNVars(scip);

   assert(sepadata->norigcuts == sepadata->nmastercuts);
   
   for( i = 0; i < sepadata->norigcuts; i++ )
   {
      /* check lhs and rhs */
      assert(SCIPisInfinity(scip, SCIProwGetRhs(sepadata->origcuts[i])) == SCIPisInfinity(scip, SCIProwGetRhs(sepadata->mastercuts[i])));
      assert(SCIPisInfinity(scip, -1 * SCIProwGetLhs(sepadata->origcuts[i])) == SCIPisInfinity(scip, -1 * SCIProwGetLhs(sepadata->mastercuts[i])));
      assert(SCIPisInfinity(scip, SCIProwGetRhs(sepadata->origcuts[i]))
         || (SCIProwGetRhs(sepadata->origcuts[i]) - SCIProwGetConstant(sepadata->origcuts[i])
            == SCIProwGetRhs(sepadata->mastercuts[i]) - SCIProwGetConstant(sepadata->mastercuts[i])));
      assert(SCIPisInfinity(scip, -1 * SCIProwGetRhs(sepadata->origcuts[i]))
         || (SCIProwGetLhs(sepadata->origcuts[i]) - SCIProwGetConstant(sepadata->origcuts[i])
            == SCIProwGetLhs(sepadata->mastercuts[i]) - SCIProwGetConstant(sepadata->mastercuts[i])));

      /* check the variables */
      norigcols = SCIProwGetNNonz(sepadata->origcuts[i]);
      origcols = SCIProwGetCols(sepadata->origcuts[i]);
      origvals = SCIProwGetVals(sepadata->origcuts[i]);
      nmastercols = SCIProwGetNNonz(sepadata->mastercuts[i]);
      mastercols = SCIProwGetCols(sepadata->mastercuts[i]);
      mastervals = SCIProwGetVals(sepadata->mastercuts[i]);

      for( v = 0; v < nmastervars; v++ )
      {
         /* get value of the variable in the master cut */
         masterval = 0.0;
         for( j = 0; j < nmastercols; j++ )
         {
            if( SCIPcolGetVar(mastercols[j]) == mastervars[v] )
            {
               masterval = mastervals[j];
               break;
            }
         }
         
         vardata = SCIPvarGetData(mastervars[v]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_MASTER);
         
         /* sum up values of corresponding original variables and compare it to the mastervalue */
         sum = 0.0;
         for( k = 0; k < vardata->data.mastervardata.norigvars; k++ )
         {
            for( l = 0; l < norigcols; l++ )
            {
               if( SCIPcolGetVar(origcols[l]) == vardata->data.mastervardata.origvars[k] )
               {
                  sum += origvals[l] * vardata->data.mastervardata.origvals[k];
               }
            }
         }
         assert(SCIPisFeasEQ(scip, masterval, sum));
         
      }
   }
   
   return SCIP_OKAY;
}
#endif




/*
 * Callback methods of separator
 */

#define sepaCopyMaster NULL

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


/** initialization method of separator (called after problem was transformed) */
#if 0
static
SCIP_DECL_SEPAINIT(sepaInitMaster)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of master separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitMaster NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 1
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
#else
#define sepaExitMaster NULL
#endif



/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_SEPAINITSOL(sepaInitsolMaster)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of master separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitsolMaster NULL
#endif


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolMaster)
{  
   SCIP* origscip;
   SCIP_SEPADATA* sepadata;
   int i;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->nmastercuts == sepadata->norigcuts);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

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
   int sum;
   SCIP_Bool feasible;

   assert(scip != NULL);
   assert(result != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPdebugMessage("sepaExeclpMaster\n");
   //SCIPdebugMessage("%d cuts are in the original LP!\n", SCIPgetNCutsApplied(origscip));
   //SCIPdebugMessage("%d cuts are in the master LP!\n", SCIPgetNCutsApplied(scip));

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
   sum = 0;
   
   cuts = SCIPgetCuts(origscip);
   ncuts = SCIPgetNCuts(origscip);

   /* save cuts in the origcuts array in the separator data */
   assert(sepadata->norigcuts == sepadata->nmastercuts);
   SCIP_CALL( ensureSizeCuts(scip, sepadata, sepadata->norigcuts + ncuts) );
   for( i = 0; i < ncuts; i++ )
   {
      sepadata->origcuts[sepadata->norigcuts] = cuts[i];
      SCIPcaptureRow(origscip, sepadata->origcuts[sepadata->norigcuts]);
      sepadata->norigcuts++;
   }


   SCIP_CALL( SCIPclearCuts(origscip) );

   mastervars = SCIPgetVars(scip);
   nmastervars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &mastervals, nmastervars) );
      
   for( i = 0; i < ncuts; i++ )
   {
      origcut = sepadata->origcuts[sepadata->norigcuts-ncuts+i];
      /* add orig cut to the original scip */
      //SCIP_CALL( SCIPaddCut(origscip, GCGrelaxGetCurrentOrigSol(origscip), origcut, FALSE) );

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
      //SCIPcaptureRow(scip, sepadata->mastercuts[sepadata->nmastercuts]);
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

#if 0
   SCIP_CALL( checkCutConsistency(scip) );
#endif

   return SCIP_OKAY;
}



/** arbitrary primal solution separation method of separator */
#if 0
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolMaster)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of master separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExecsolMaster NULL
#endif




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

   SCIP_CALL( SCIPallocMemoryArray(scip, &(sepadata->origcuts), STARTMAXCUTS) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sepadata->mastercuts), STARTMAXCUTS) );
   sepadata->maxcuts = STARTMAXCUTS;
   sepadata->norigcuts = 0;
   sepadata->nmastercuts = 0;


   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaCopyMaster, sepaFreeMaster, sepaInitMaster, sepaExitMaster, 
         sepaInitsolMaster, sepaExitsolMaster,
         sepaExeclpMaster, sepaExecsolMaster,
         sepadata) );

   /* add master separator parameters */
   /* TODO: (optional) add separator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}


/* returns the array of original cuts saved in the separator data */
SCIP_ROW** GCGsepaGetOrigcuts(
   SCIP*                 scip
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
   SCIP*                 scip
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
   SCIP*                 scip
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
   SCIP*                 scip
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
