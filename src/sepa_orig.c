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

/**@file   sepa_orig.c
 * @ingroup SEPARATORS
 * @brief  orig separator
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "sepa_orig.h"
#include "probdata_gcg.h"


#define SEPA_NAME              "orig"
#define SEPA_DESC              "separator for gcg separating cuts in the original problem"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                     1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */




/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
};




/*
 * Local methods
 */


/*
 * main separation method
 */

/** searches and adds clique cuts that separate the given primal solution */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< the cut separator itself */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{  
   SCIP* origprob;
   SCIP_VARDATA* vardata;
   SCIP_SOL* origsol;
   SCIP_VAR** mastervars;
   int nmastervars;
   SCIP_Bool delayed;
   SCIP_Bool cutoff;
   int v;
   int i;

   assert(scip != NULL);

   origprob = GCGprobGetOrigprob(scip);
   assert(origprob != NULL);

   *result = SCIP_DIDNOTRUN;

   if ( sol == NULL )
   {
      SCIP_CALL( SCIPcreateLPSol(scip, &sol, NULL) );
   }
   //SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );

   nmastervars = SCIPgetNVars(scip);
   mastervars = SCIPgetVars(scip);

   /* create solution for the original problem and conververt current solution to this solution space */
   SCIP_CALL( SCIPcreateSol(origprob, &origsol, NULL) );

   for ( v = 0; v < nmastervars; v++ )
   {
      if ( !SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, mastervars[v])) )
      {
         vardata = SCIPvarGetData(mastervars[v]);
         assert(vardata != NULL);
         assert(vardata->vartype == GCG_VARTYPE_MASTER);
         assert(vardata->data.mastervardata.norigvars > 0);
         for ( i = 0; i < vardata->data.mastervardata.norigvars; i++ )
         {
            SCIP_CALL( SCIPincSolVal(origprob, origsol, vardata->data.mastervardata.origvars[i],
                  vardata->data.mastervardata.origvals[i]*SCIPgetSolVal(scip, sol, mastervars[v])) );
         }
      }
   }
   
   //SCIP_CALL( SCIPprintSol(origprob, origsol, NULL, FALSE) );

   assert(SCIPgetNCuts(origprob) == 0);

   SCIP_CALL( SCIPseparateSol(origprob, origsol, TRUE, FALSE, &delayed, &cutoff) );

   printf("%d cuts found!\n", SCIPgetNCuts(origprob));

   SCIP_CALL( SCIPprintStatistics(origprob, NULL) );

   return SCIP_OKAY;
}





/*
 * Callback methods of separator
 */

/** destructor of separator to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_SEPAFREE(sepaFreeOrig)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of orig separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaFreeOrig NULL
#endif


/** initialization method of separator (called after problem was transformed) */
#if 0
static
SCIP_DECL_SEPAINIT(sepaInitOrig)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of orig separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitOrig NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_SEPAEXIT(sepaExitOrig)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of orig separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitOrig NULL
#endif


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_SEPAINITSOL(sepaInitsolOrig)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of orig separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitsolOrig NULL
#endif


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolOrig)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of orig separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitsolOrig NULL
#endif


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpOrig)
{  
   printf("sepaexeclporig\n");
   /* separate cuts on the LP solution */
   SCIP_CALL( separateCuts(scip, sepa, NULL, result) );

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolOrig)
{  
   printf("sepaexecsolorig\n");
   /* separate cuts on the given primal solution */
   SCIP_CALL( separateCuts(scip, sepa, sol, result) );

   return SCIP_OKAY;
}



/*
 * separator specific interface methods
 */

/** creates the orig separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaOrig(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create orig separator data */
   sepadata = NULL;

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_DELAY,
         sepaFreeOrig, sepaInitOrig, sepaExitOrig, 
         sepaInitsolOrig, sepaExitsolOrig,
         sepaExeclpOrig, sepaExecsolOrig,
         sepadata) );

   /* add orig separator parameters */

   return SCIP_OKAY;
}
