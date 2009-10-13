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

/**@file   disp_gcg.c
 * @ingroup DISPLAYS
 * @brief  gcg display columns
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "disp_gcg.h"


#define DISP_NAME_MLPITERATIONS  "mlpiterations"
#define DISP_DESC_MLPITERATIONS  "number of simplex iterations in the master"
#define DISP_HEAD_MLPITERATIONS  "MLP iter"
#define DISP_WIDT_MLPITERATIONS  8
#define DISP_PRIO_MLPITERATIONS  100000
#define DISP_POSI_MLPITERATIONS  8000
#define DISP_STRI_MLPITERATIONS  TRUE

#define DISP_NAME_MEMUSED       "totalmemused"
#define DISP_DESC_MEMUSED       "total number of bytes used in block memory"
#define DISP_HEAD_MEMUSED       "tmem"
#define DISP_WIDT_MEMUSED       5
#define DISP_PRIO_MEMUSED       1000000
#define DISP_POSI_MEMUSED       1500
#define DISP_STRI_MEMUSED       TRUE

#define DISP_NAME_MVARS         "mvars"
#define DISP_DESC_MVARS         "number of variables in the master problem"
#define DISP_HEAD_MVARS         "mvars"
#define DISP_WIDT_MVARS         5
#define DISP_PRIO_MVARS         100000
#define DISP_POSI_MVARS         8100
#define DISP_STRI_MVARS         TRUE

#define DISP_NAME_MCONSS        "mconss"
#define DISP_DESC_MCONSS        "number of globally valid constraints in the master problem"
#define DISP_HEAD_MCONSS        "mcons"
#define DISP_WIDT_MCONSS        5
#define DISP_PRIO_MCONSS        100000
#define DISP_POSI_MCONSS        8200
#define DISP_STRI_MCONSS        TRUE

#define DISP_NAME_MCUTS         "mcuts"
#define DISP_DESC_MCUTS         "total number of cuts applied to the master LPs"
#define DISP_HEAD_MCUTS         "mcuts"
#define DISP_WIDT_MCUTS         5
#define DISP_PRIO_MCUTS         100000
#define DISP_POSI_MCUTS         8400
#define DISP_STRI_MCUTS         TRUE


/*
 * Callback methods
 */

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputMlpiterations)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MLPITERATIONS) == 0);
   assert(scip != NULL);

   SCIPdispLongint(file, SCIPgetNLPIterations(GCGrelaxGetMasterprob(scip)), DISP_WIDT_MLPITERATIONS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputMemused)
{  /*lint --e{715}*/
   SCIP_Longint memused;
   int i;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MEMUSED) == 0);
   assert(scip != NULL);

   memused = SCIPgetMemUsed(scip);
   memused += SCIPgetMemUsed(GCGrelaxGetMasterprob(scip));
   for ( i = 0; i < GCGrelaxGetNPricingprobs(scip); i++ )
   {
      memused += SCIPgetMemUsed(GCGrelaxGetPricingprob(scip, i));
   }

   SCIPdispLongint(file, memused, DISP_WIDT_MEMUSED);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputMvars)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MVARS) == 0);
   assert(scip != NULL);

   SCIPdispInt(file, SCIPgetNVars(GCGrelaxGetMasterprob(scip)), DISP_WIDT_MVARS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputMconss)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MCONSS) == 0);
   assert(scip != NULL);

   SCIPdispInt(file, SCIPgetNConss(GCGrelaxGetMasterprob(scip)), DISP_WIDT_MCONSS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputMcuts)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MCUTS) == 0);
   assert(scip != NULL);

   SCIPdispInt(file, SCIPgetNCutsApplied(GCGrelaxGetMasterprob(scip)), DISP_WIDT_MCUTS);

   return SCIP_OKAY;
}

/*
 * gcg display columns specific interface methods
 */

/** includes the gcg display columns in SCIP */
SCIP_RETCODE SCIPincludeDispGcg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_MLPITERATIONS, DISP_DESC_MLPITERATIONS, DISP_HEAD_MLPITERATIONS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMlpiterations, NULL, 
         DISP_WIDT_MLPITERATIONS, DISP_PRIO_MLPITERATIONS, DISP_POSI_MLPITERATIONS, DISP_STRI_MLPITERATIONS) );
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_MEMUSED, DISP_DESC_MEMUSED, DISP_HEAD_MEMUSED,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMemused, NULL, 
         DISP_WIDT_MEMUSED, DISP_PRIO_MEMUSED, DISP_POSI_MEMUSED, DISP_STRI_MEMUSED) );
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_MVARS, DISP_DESC_MVARS, DISP_HEAD_MVARS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMvars, NULL, 
         DISP_WIDT_MVARS, DISP_PRIO_MVARS, DISP_POSI_MVARS, DISP_STRI_MVARS) );
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_MCONSS, DISP_DESC_MCONSS, DISP_HEAD_MCONSS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMconss, NULL, 
         DISP_WIDT_MCONSS, DISP_PRIO_MCONSS, DISP_POSI_MCONSS, DISP_STRI_MCONSS) );
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_MCUTS, DISP_DESC_MCUTS, DISP_HEAD_MCUTS,
         SCIP_DISPSTATUS_AUTO, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMcuts, NULL, 
         DISP_WIDT_MCUTS, DISP_PRIO_MCUTS, DISP_POSI_MCUTS, DISP_STRI_MCUTS) );

   return SCIP_OKAY;
}

