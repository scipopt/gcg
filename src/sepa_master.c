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
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */


#define STARTMAXCUTS 50
#define MAXCUTSINC 20


/*
 * Data structures
 */

/* TODO: fill in the necessary separator data */

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

   if ( sepadata->maxcuts < size )
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
   SCIP_ROW* mastercut;
   SCIP_ROW* origcut;
   SCIP_COL** cols;
   SCIP_VAR** rowvars;
   int ncols;
   SCIP_Real* vals;
   SCIP_SEPADATA* sepadata;
   FILE* origcutsfile;

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

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPdebugMessage("sepaExeclpMaster\n");
   SCIPdebugMessage("%d cuts are in the original LP!\n", SCIPgetNCutsApplied(origscip));
   SCIPdebugMessage("%d cuts are in the master LP!\n", SCIPgetNCutsApplied(scip));

   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( GCGrelaxUpdateCurrentSol(origscip) );

   SCIP_CALL( SCIPseparateSol(origscip, GCGrelaxGetCurrentOrigSol(origscip),
         FALSE, FALSE, &delayed, &cutoff) );   

   SCIPdebugMessage("SCIPseparateSol() found %d cuts!\n", SCIPgetNCuts(origscip));
   sum = 0;
   
   cuts = SCIPgetCuts(origscip);
   ncuts = SCIPgetNCuts(origscip);

   origcutsfile = fopen("origcuts.txt", "w");
   //printf("origcutsfile %s NULL\n", (origcutsfile == NULL ? "==" : "!="));
   

   /* save cuts in the origcuts array in the separator data */
   assert(sepadata->norigcuts == sepadata->nmastercuts);
   SCIP_CALL( ensureSizeCuts(scip, sepadata, sepadata->norigcuts + ncuts) );
   for ( i = 0; i < ncuts; i++ )
   {
      sepadata->origcuts[sepadata->norigcuts] = cuts[i];
      SCIProwCapture(sepadata->origcuts[sepadata->norigcuts]);
      sepadata->norigcuts++;
   }


   SCIP_CALL( SCIPclearCuts(origscip) );

   mastervars = SCIPgetVars(scip);
   nmastervars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &mastervals, nmastervars) );
      
   for ( i = 0; i < ncuts; i++ )
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
      for ( j = 0; j < ncols; j++ )
      {
         rowvars[j] = SCIPcolGetVar(cols[j]);
      }

      /* create new cut in the master problem */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "m_%s", SCIProwGetName(origcut));
      SCIP_CALL( SCIPcreateEmptyRow(scip, &mastercut, name, 
            ( SCIPisInfinity(scip, -SCIProwGetLhs(origcut)) ? 
               SCIProwGetLhs(origcut) : SCIProwGetLhs(origcut) - SCIProwGetConstant(origcut)), 
            ( SCIPisInfinity(scip, SCIProwGetRhs(origcut)) ? 
               SCIProwGetRhs(origcut) : SCIProwGetRhs(origcut) - SCIProwGetConstant(origcut)), 
            SCIProwIsLocal(origcut), TRUE, FALSE) );

      /* transform the original variables to master variables and add them to the cut */
      GCGrelaxTransformOrigvalsToMastervals(scip, rowvars, vals, ncols, mastervars, mastervals, nmastervars);
      SCIP_CALL( SCIPaddVarsToRow(scip, mastercut, nmastervars, mastervars, mastervals) );

      /* add the cut to the master problem */
      SCIP_CALL( SCIPaddCut(scip, NULL, mastercut, FALSE) );
      sepadata->mastercuts[sepadata->nmastercuts] = mastercut;
      SCIProwCapture(sepadata->mastercuts[sepadata->nmastercuts]);
      sepadata->nmastercuts++;

#ifdef SCIP_DEBUG
      //SCIPdebugMessage("Cut %d:\n", i);
      SCIP_CALL( SCIPprintRow(origscip, origcut, origcutsfile) );
      //SCIP_CALL( SCIPprintRow(scip, mastercut, NULL) );
      //SCIPdebugMessage("\n\n");
#endif 

      SCIPfreeBufferArray(scip, &rowvars);
      
   }
   
   fclose(origcutsfile);

   SCIPdebugMessage("%d cuts are in the original sepastore!\n", SCIPgetNCuts(origscip));
   SCIPdebugMessage("%d cuts are in the master sepastore!\n", SCIPgetNCuts(scip));

   SCIPfreeBufferArray(scip, &mastervals);

   assert(sepadata->norigcuts == sepadata->nmastercuts );

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
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_DELAY,
         sepaFreeMaster, sepaInitMaster, sepaExitMaster, 
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
