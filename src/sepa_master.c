/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"
#define SCIP_DEBUG
/**@file   sepa_master.c
 * @ingroup SEPARATORS
 * @brief  master separator
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

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
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */




/*
 * Data structures
 */

/* TODO: fill in the necessary separator data */

/** separator data */
struct SCIP_SepaData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of separator
 */

/* TODO: Implement all necessary separator methods. The methods with an #if 0 ... #else #define ... are optional */

/** destructor of separator to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_SEPAFREE(sepaFreeMaster)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of master separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaFreeMaster NULL
#endif


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
#if 0
static
SCIP_DECL_SEPAEXIT(sepaExitMaster)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of master separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

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
#if 0
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolMaster)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of master separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitsolMaster NULL
#endif


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpMaster)
{  
   SCIP*   origscip;
   SCIP_Bool delayed;
   SCIP_Bool cutoff;
   SCIP_ROW** cuts;
   SCIP_ROW** cutscopy;
   SCIP_ROW* mastercut;
   SCIP_COL** cols;
   SCIP_VAR** rowvars;
   int ncols;
   SCIP_Real* vals;

   SCIP_VAR** mastervars;
   SCIP_Real* mastervals;
   int nmastervars;

   char name[SCIP_MAXSTRLEN];

   int ncuts;
   int i;
   int j;
   int sum;

   assert(scip != NULL);
   assert(result != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   SCIPdebugMessage("sepaExeclpMaster\n");

   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( GCGrelaxUpdateCurrentSol(origscip) );

   SCIP_CALL( SCIPseparateSol(origscip, GCGrelaxGetCurrentOrigSol(origscip),
         FALSE, FALSE, &delayed, &cutoff) );   

   SCIPdebugMessage("SCIPseparateSol() found %d cuts!\n", SCIPgetNCuts(origscip));
   sum = 0;
   
   cuts = SCIPgetCuts(origscip);
   ncuts = SCIPgetNCuts(origscip);
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &cutscopy, cuts, ncuts) );
   for ( i = 0; i < ncuts; i++ )
   {
      SCIProwCapture(cutscopy[i]);
   }



   SCIP_CALL( SCIPclearCuts(origscip) );

   mastervars = SCIPgetVars(scip);
   nmastervars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &mastervals, nmastervars) );
   
   
   for ( i = 0; i < ncuts; i++ )
   {
      SCIP_CALL( SCIPaddCut(origscip, GCGrelaxGetCurrentOrigSol(origscip), cutscopy[i], TRUE) );

      ncols = SCIProwGetNNonz(cutscopy[i]);
      cols = SCIProwGetCols(cutscopy[i]);
      vals = SCIProwGetVals(cutscopy[i]);

      SCIP_CALL( SCIPallocBufferArray(scip, &rowvars, ncols) );
      for ( j = 0; j < ncols; j++ )
      {
         rowvars[j] = SCIPcolGetVar(cols[j]);
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "m_%s", SCIProwGetName(cutscopy[i]));

      assert(SCIPisFeasZero(scip, SCIProwGetConstant(cutscopy[i])));
      SCIP_CALL( SCIPcreateEmptyRow(scip, &mastercut, name, SCIProwGetLhs(cutscopy[i]), SCIProwGetRhs(cutscopy[i]), 
            SCIProwIsLocal(cutscopy[i]), TRUE, FALSE) );

      GCGrelaxTransformOrigvalsToMastervals(rowvars, vals, ncols, mastervars, mastervals, nmastervars);

      for ( j = 0; j < nmastervars; j++ )
      {
         SCIP_CALL( SCIPaddVarsToRow(scip, mastercut, nmastervars, mastervars, mastervals) );
      }

      SCIP_CALL( SCIPaddCut(scip, NULL, mastercut, TRUE) );

      SCIPdebugMessage("Cut %d:\n", i);
      SCIP_CALL( SCIPprintRow(origscip, cutscopy[i], NULL) );
      SCIP_CALL( SCIPprintRow(scip, mastercut, NULL) );
      SCIPdebugMessage("\n\n");

      SCIPfreeBufferArray(scip, &rowvars);
      
   }

   SCIPdebugMessage("%d cuts are in the sepastore!\n", SCIPgetNCuts(origscip));
   SCIPdebugMessage("%d cuts are in the LP!\n", SCIPgetNCutsApplied(origscip));

   SCIPfreeBufferArray(scip, &mastervals);
   SCIPfreeMemoryArray(scip, &cutscopy);

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
   sepadata = NULL;
   /* TODO: (optional) create separator specific data here */

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_DELAY,
         sepaFreeMaster, sepaInitMaster, sepaExitMaster, 
         sepaInitsolMaster, sepaExitsolMaster,
         sepaExeclpMaster, sepaExecsolMaster,
         sepadata) );

   /* add master separator parameters */
   /* TODO: (optional) add separator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
