/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2012 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   disp_gcg.c
 * @ingroup DISPLAYS
 * @brief  GCG display columns
 * @author Gerald Gamrath
 * @author Christian Puchert
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "disp_gcg.h"
#include "relax_gcg.h"

#define DISP_NAME_SOLFOUND      "solfound"
#define DISP_DESC_SOLFOUND      "letter that indicates the heuristic, that found the solution"
#define DISP_HEAD_SOLFOUND      "  "
#define DISP_WIDT_SOLFOUND      2
#define DISP_PRIO_SOLFOUND      80000
#define DISP_POSI_SOLFOUND      0
#define DISP_STRI_SOLFOUND      FALSE

#define DISP_NAME_MLPITERATIONS  "mlpiterations"
#define DISP_DESC_MLPITERATIONS  "number of simplex iterations in the master"
#define DISP_HEAD_MLPITERATIONS  "MLP iter"
#define DISP_WIDT_MLPITERATIONS  8
#define DISP_PRIO_MLPITERATIONS  80000
#define DISP_POSI_MLPITERATIONS  1001
#define DISP_STRI_MLPITERATIONS  TRUE

#define DISP_NAME_LPAVGITERS    "lpavgiterations"
#define DISP_DESC_LPAVGITERS    "average number of LP iterations since the last output line"
#define DISP_HEAD_LPAVGITERS    "LP it/n"
#define DISP_WIDT_LPAVGITERS    7
#define DISP_PRIO_LPAVGITERS    500
#define DISP_POSI_LPAVGITERS    1400
#define DISP_STRI_LPAVGITERS    TRUE

#define DISP_NAME_MLPAVGITERS    "mlpavgiterations"
#define DISP_DESC_MLPAVGITERS    "average number of LP iterations in the master"
#define DISP_HEAD_MLPAVGITERS    "MLP it/n"
#define DISP_WIDT_MLPAVGITERS    8
#define DISP_PRIO_MLPAVGITERS    25000
#define DISP_POSI_MLPAVGITERS    1401
#define DISP_STRI_MLPAVGITERS    TRUE

#define DISP_NAME_MLPCOND        "mlpcond"
#define DISP_DESC_MLPCOND        "estimate on condition number of LP master solution"
#define DISP_HEAD_MLPCOND        "MLP cond"
#define DISP_WIDT_MLPCOND        8
#define DISP_PRIO_MLPCOND        0
#define DISP_POSI_MLPCOND        1451
#define DISP_STRI_MLPCOND        TRUE

#define DISP_NAME_MEMUSED       "memused"
#define DISP_DESC_MEMUSED       "total number of bytes used in block memory"
#define DISP_HEAD_MEMUSED       "mem"
#define DISP_WIDT_MEMUSED       5
#define DISP_PRIO_MEMUSED       20000
#define DISP_POSI_MEMUSED       1500
#define DISP_STRI_MEMUSED       TRUE

#define DISP_NAME_VARS          "vars"
#define DISP_DESC_VARS          "number of variables in the original problem"
#define DISP_HEAD_VARS          "ovars"
#define DISP_WIDT_VARS          5
#define DISP_PRIO_VARS          3000
#define DISP_POSI_VARS          3000
#define DISP_STRI_VARS          TRUE

#define DISP_NAME_CONSS         "conss"
#define DISP_DESC_CONSS         "number of globally valid constraints in the problem"
#define DISP_HEAD_CONSS         "ocons"
#define DISP_WIDT_CONSS         5
#define DISP_PRIO_CONSS         3100
#define DISP_POSI_CONSS         3100
#define DISP_STRI_CONSS         TRUE

#define DISP_NAME_CUTS          "cuts"
#define DISP_DESC_CUTS          "total number of cuts applied to the original LPs"
#define DISP_HEAD_CUTS          "ocuts"
#define DISP_WIDT_CUTS          5
#define DISP_PRIO_CUTS          100
#define DISP_POSI_CUTS          3500
#define DISP_STRI_CUTS          TRUE

#define DISP_NAME_SEPAROUNDS    "separounds"
#define DISP_DESC_SEPAROUNDS    "number of separation rounds performed at the current node"
#define DISP_HEAD_SEPAROUNDS    "sepa"
#define DISP_WIDT_SEPAROUNDS    4
#define DISP_PRIO_SEPAROUNDS    100
#define DISP_POSI_SEPAROUNDS    3600
#define DISP_STRI_SEPAROUNDS    TRUE

#define DISP_NAME_POOLSIZE      "poolsize"
#define DISP_DESC_POOLSIZE      "number of LP rows in the cut pool"
#define DISP_HEAD_POOLSIZE      "pool"
#define DISP_WIDT_POOLSIZE      5
#define DISP_PRIO_POOLSIZE      50
#define DISP_POSI_POOLSIZE      3700
#define DISP_STRI_POOLSIZE      TRUE

#define DISP_NAME_LPOBJ         "lpobj"
#define DISP_DESC_LPOBJ         "current LP objective value"
#define DISP_HEAD_LPOBJ         "lpobj"
#define DISP_WIDT_LPOBJ         14
#define DISP_PRIO_LPOBJ         300
#define DISP_POSI_LPOBJ         6500
#define DISP_STRI_LPOBJ         TRUE

#define DISP_NAME_MVARS         "mvars"
#define DISP_DESC_MVARS         "number of variables in the master problem"
#define DISP_HEAD_MVARS         "mvars"
#define DISP_WIDT_MVARS         5
#define DISP_PRIO_MVARS         70000
#define DISP_POSI_MVARS         3050
#define DISP_STRI_MVARS         TRUE

#define DISP_NAME_MCONSS        "mconss"
#define DISP_DESC_MCONSS        "number of globally valid constraints in the master problem"
#define DISP_HEAD_MCONSS        "mcons"
#define DISP_WIDT_MCONSS        5
#define DISP_PRIO_MCONSS        70000
#define DISP_POSI_MCONSS        3150
#define DISP_STRI_MCONSS        TRUE

#define DISP_NAME_MCUTS         "mcuts"
#define DISP_DESC_MCUTS         "total number of cuts applied to the master LPs"
#define DISP_HEAD_MCUTS         "mcuts"
#define DISP_WIDT_MCUTS         5
#define DISP_PRIO_MCUTS         80000
#define DISP_POSI_MCUTS         3550
#define DISP_STRI_MCUTS         TRUE

/**TODO:
 *
 * Degeneracy: in %
 */

/*
 * Callback methods
 */

/** copy method for display plugins (called when SCIP copies plugins) */
static
SCIP_DECL_DISPCOPY(dispCopyDefault)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(disp != NULL);

   /* call inclusion method of dialog */
   SCIP_CALL( SCIPincludeDispGcg(scip) );

   return SCIP_OKAY;
}

/** solving process initialization method of display column (called when branch and bound process is about to begin) */
static
SCIP_DECL_DISPINITSOL(SCIPdispInitsolSolFound)
{  /*lint --e{715}*/

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_SOLFOUND) == 0);
   assert(scip != NULL);

   SCIPdispSetData(disp, (SCIP_DISPDATA*)SCIPgetBestSol(scip));

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' for character of best solution */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputSolFound)
{  /*lint --e{715}*/
   SCIP* masterprob;
   SCIP_SOL* origsol;
   SCIP_SOL* mastersol;
   SCIP_DISPDATA* dispdata;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_SOLFOUND) == 0);
   assert(scip != NULL);

   /* get master problem */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   origsol = SCIPgetBestSol(scip);
   if( origsol == NULL )
      SCIPdispSetData(disp, NULL);

   if( SCIPgetStage(masterprob) >= SCIP_STAGE_SOLVING )
      mastersol = SCIPgetBestSol(masterprob);
   else
      mastersol = NULL;

   dispdata = SCIPdispGetData(disp);
   if( origsol != (SCIP_SOL*)dispdata )
   {
      SCIPinfoMessage(scip, file, "%c", (SCIPgetSolHeur(scip, origsol) == NULL ? '*'
            : SCIPheurGetDispchar(SCIPgetSolHeur(scip, origsol))));
      /* If the solution was obtained in the master problem, display whether it came from its
       * LP relaxation or from the master heuristics */
      if( SCIPgetSolHeur(scip, origsol) == NULL && (mastersol != NULL) )
      {
         SCIPinfoMessage(scip, file, "%c", (SCIPgetSolHeur(masterprob, mastersol) == NULL ? '*'
               : SCIPheurGetDispchar(SCIPgetSolHeur(masterprob, mastersol))));
      }
      else
      {
         SCIPinfoMessage(scip, file, " ");
      }
      SCIPdispSetData(disp, (SCIP_DISPDATA*)origsol);
   }
   else
      SCIPinfoMessage(scip, file, "  ");

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' for number of master LP iterations */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputMlpiterations)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MLPITERATIONS) == 0);
   assert(scip != NULL);

   if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) >= SCIP_STAGE_SOLVING )
   {
      SCIPdispLongint(SCIPgetMessagehdlr(scip), file, SCIPgetNLPIterations(GCGrelaxGetMasterprob(scip)), DISP_WIDT_MLPITERATIONS);
   }
   else
   {
      SCIPdispLongint(SCIPgetMessagehdlr(scip), file, 0LL, DISP_WIDT_MLPITERATIONS);
   }

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' for number of average LP iterations */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputNLPAvgIters)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_LPAVGITERS) == 0);
   assert(scip != NULL);

   /**@todo Currently we are using the total number of nodes to compute the average LP iterations number. The reason for
    *       that is, that for the LP iterations only the total number (over all runs) are stored in the statistics. It
    *       would be nicer if the statistic also stores the number of LP iterations for the current run similar to the
    *       nodes.
    */

   if( SCIPgetNNodes(scip) < 2 )
      SCIPinfoMessage(scip, file, "     - ");
   else
      SCIPinfoMessage(scip, file, "%6.1f ",
         (SCIPgetNLPIterations(scip) - SCIPgetNRootLPIterations(scip)) / (SCIP_Real)(SCIPgetNTotalNodes(scip) - 1) );

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' for number of average master LP iterations */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputNMLPAvgIters)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MLPAVGITERS) == 0);
   assert(scip != NULL);

   /**@todo Currently we are using the total number of nodes to compute the average LP iterations number. The reason for
    *       that is, that for the LP iterations only the total number (over all runs) are stored in the statistics. It
    *       would be nicer if the statistic also stores the number of LP iterations for the current run similar to the
    *       nodes.
    */

   if( SCIPgetNNodes(scip) < 2 )
      SCIPinfoMessage(scip, file, "     - ");
   else
      SCIPinfoMessage(scip, file, "%6.1f ",
         (SCIPgetNLPIterations(GCGrelaxGetMasterprob(scip)) - SCIPgetNRootLPIterations(GCGrelaxGetMasterprob(scip)))
         / (SCIP_Real)(SCIPgetNNodes(GCGrelaxGetMasterprob(scip)) - 1) );

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' for estimate on master LP condition */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputMLPCondition)
{  /*lint --e{715}*/
   SCIP_LPI* lpi;
   SCIP_Real cond;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_LPCOND) == 0);
   assert(scip != NULL);

   SCIP_CALL( SCIPgetLPI(GCGrelaxGetMasterprob(scip), &lpi) );
   if( lpi == NULL )
   {
      SCIPinfoMessage(scip, file, "     - ");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPlpiGetRealSolQuality(lpi, SCIP_LPSOLQUALITY_ESTIMCONDITION, &cond) );

   if( cond == SCIP_INVALID )  /*lint !e777*/
      SCIPinfoMessage(scip, file, "   n/a ", cond);
   else
      SCIPinfoMessage(scip, file, "%.1e", cond);

   return SCIP_OKAY;
}


/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputMemUsed)
{  /*lint --e{715}*/
   SCIP_Longint memused;
   int i;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MEMUSED) == 0);
   assert(scip != NULL);

   memused = SCIPgetMemUsed(scip);
   memused += SCIPgetMemUsed(GCGrelaxGetMasterprob(scip));
   for( i = 0; i < GCGrelaxGetNPricingprobs(scip); i++ )
   {
      memused += SCIPgetMemUsed(GCGrelaxGetPricingprob(scip, i));
   }

   SCIPdispLongint(SCIPgetMessagehdlr(scip), file, memused, DISP_WIDT_MEMUSED);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' for number of variables */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputNVars)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_VARS) == 0);
   assert(scip != NULL);

   SCIPdispInt(SCIPgetMessagehdlr(scip), file, SCIPgetNVars(scip), DISP_WIDT_VARS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' for number of constraints */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputNConss)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_CONSS) == 0);
   assert(scip != NULL);

   SCIPdispInt(SCIPgetMessagehdlr(scip), file, SCIPgetNConss(scip), DISP_WIDT_CONSS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' for number of applied cuts */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputNAppliedCuts)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_CUTS) == 0);
   assert(scip != NULL);

   SCIPdispInt(SCIPgetMessagehdlr(scip), file, SCIPgetNCutsApplied(scip), DISP_WIDT_CUTS);

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' for number of separation rounds */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputNSepaRounds)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_SEPAROUNDS) == 0);
   assert(scip != NULL);

   if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) == SCIP_STAGE_SOLVING )
   {
      SCIPdispInt(SCIPgetMessagehdlr(scip), file, SCIPgetNSepaRounds(GCGrelaxGetMasterprob(scip)), DISP_WIDT_SEPAROUNDS);
   }
   else
   {
      SCIPdispInt(SCIPgetMessagehdlr(scip), file, 0, DISP_WIDT_SEPAROUNDS);
   }

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' for number of current rows in the cut pool */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputCutPoolSize)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_POOLSIZE) == 0);
   assert(scip != NULL);

   if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) >= SCIP_STAGE_SOLVING )
   {
      SCIPdispInt(SCIPgetMessagehdlr(scip), file, SCIPgetNPoolCuts(GCGrelaxGetMasterprob(scip)), DISP_WIDT_POOLSIZE);
   }
   else
   {
      SCIPdispInt(SCIPgetMessagehdlr(scip), file, 0, DISP_WIDT_POOLSIZE);
   }

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' for LP objective value */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputLPObjval)
{  /*lint --e{715}*/
   SCIP_Real lpobj;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_LPOBJ) == 0);
   assert(scip != NULL);

   if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) != SCIP_STAGE_SOLVING || SCIPgetLPSolstat(GCGrelaxGetMasterprob(scip)) == SCIP_LPSOLSTAT_NOTSOLVED )
   {
      SCIPinfoMessage(scip, file, "      --      ");
   }
   else
   {
      lpobj = SCIPgetLPObjval(GCGrelaxGetMasterprob(scip));
      if( SCIPisInfinity(scip, -lpobj) )
         SCIPinfoMessage(scip, file, "      --      ");
      else if( SCIPisInfinity(scip, lpobj) )
         SCIPinfoMessage(scip, file, "    cutoff    ");
      else
         SCIPinfoMessage(scip, file, "%13.6e ", lpobj);
   }

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputMvars)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MVARS) == 0);
   assert(scip != NULL);

   if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) >= SCIP_STAGE_SOLVING )
   {
      SCIPdispInt(SCIPgetMessagehdlr(scip), file, SCIPgetNVars(GCGrelaxGetMasterprob(scip)), DISP_WIDT_MVARS);
   }
   else
   {
      SCIPdispInt(SCIPgetMessagehdlr(scip), file, 0, DISP_WIDT_MVARS);
   }

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputMconss)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MCONSS) == 0);
   assert(scip != NULL);

   if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) >= SCIP_STAGE_SOLVING )
   {
      SCIPdispInt(SCIPgetMessagehdlr(scip), file, SCIPgetNConss(GCGrelaxGetMasterprob(scip)), DISP_WIDT_MCONSS);
   }
   else
   {
      SCIPdispInt(SCIPgetMessagehdlr(scip), file, 0, DISP_WIDT_MCONSS);
   }

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(SCIPdispOutputMcuts)
{  /*lint --e{715}*/
   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_NAME_MCUTS) == 0);
   assert(scip != NULL);

   if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) >= SCIP_STAGE_SOLVING )
   {
      SCIPdispInt(SCIPgetMessagehdlr(scip), file, SCIPgetNCutsApplied(GCGrelaxGetMasterprob(scip)), DISP_WIDT_MCUTS);
   }
   else
   {
      SCIPdispInt(SCIPgetMessagehdlr(scip), file, 0, DISP_WIDT_MCUTS);
   }



   return SCIP_OKAY;
}

/*
 * default display columns specific interface methods
 */

/** includes the default display columns in SCIP */
SCIP_RETCODE SCIPincludeDispGcg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DISP* tmpdisp;

   tmpdisp = SCIPfindDisp(scip, DISP_NAME_SOLFOUND);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_SOLFOUND, DISP_DESC_SOLFOUND, DISP_HEAD_SOLFOUND,
            SCIP_DISPSTATUS_AUTO,
            dispCopyDefault,
            NULL, NULL, NULL, SCIPdispInitsolSolFound, NULL, SCIPdispOutputSolFound, NULL,
            DISP_WIDT_SOLFOUND, DISP_PRIO_SOLFOUND, DISP_POSI_SOLFOUND, DISP_STRI_SOLFOUND) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_MLPITERATIONS);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_MLPITERATIONS, DISP_DESC_MLPITERATIONS, DISP_HEAD_MLPITERATIONS,
            SCIP_DISPSTATUS_AUTO,
            dispCopyDefault,
            NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMlpiterations, NULL,
            DISP_WIDT_MLPITERATIONS, DISP_PRIO_MLPITERATIONS, DISP_POSI_MLPITERATIONS, DISP_STRI_MLPITERATIONS) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_MLPAVGITERS);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_MLPAVGITERS, DISP_DESC_MLPAVGITERS, DISP_HEAD_MLPAVGITERS,
            SCIP_DISPSTATUS_AUTO,
            dispCopyDefault,
            NULL, NULL, NULL, NULL, NULL, SCIPdispOutputNMLPAvgIters, NULL,
            DISP_WIDT_MLPAVGITERS, DISP_PRIO_MLPAVGITERS, DISP_POSI_MLPAVGITERS, DISP_STRI_MLPAVGITERS) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_LPAVGITERS);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_LPAVGITERS, DISP_DESC_LPAVGITERS, DISP_HEAD_LPAVGITERS,
            SCIP_DISPSTATUS_AUTO,
            dispCopyDefault,
            NULL, NULL, NULL, NULL, NULL, SCIPdispOutputNLPAvgIters, NULL,
            DISP_WIDT_LPAVGITERS, DISP_PRIO_LPAVGITERS, DISP_POSI_LPAVGITERS, DISP_STRI_LPAVGITERS) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_MLPCOND);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_MLPCOND, DISP_DESC_MLPCOND, DISP_HEAD_MLPCOND,
            SCIP_DISPSTATUS_AUTO,
            dispCopyDefault,
            NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMLPCondition, NULL,
            DISP_WIDT_MLPCOND, DISP_PRIO_MLPCOND, DISP_POSI_MLPCOND, DISP_STRI_MLPCOND) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_MEMUSED);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_MEMUSED, DISP_DESC_MEMUSED, DISP_HEAD_MEMUSED,
            SCIP_DISPSTATUS_AUTO,
            dispCopyDefault,
            NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMemUsed, NULL,
            DISP_WIDT_MEMUSED, DISP_PRIO_MEMUSED, DISP_POSI_MEMUSED, DISP_STRI_MEMUSED) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_VARS);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_VARS, DISP_DESC_VARS, DISP_HEAD_VARS,
            SCIP_DISPSTATUS_AUTO,
            dispCopyDefault,
            NULL, NULL, NULL, NULL, NULL, SCIPdispOutputNVars, NULL,
            DISP_WIDT_VARS, DISP_PRIO_VARS, DISP_POSI_VARS, DISP_STRI_VARS) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_CONSS);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_CONSS, DISP_DESC_CONSS, DISP_HEAD_CONSS,
            SCIP_DISPSTATUS_AUTO,
            dispCopyDefault,
            NULL, NULL, NULL, NULL, NULL, SCIPdispOutputNConss, NULL,
            DISP_WIDT_CONSS, DISP_PRIO_CONSS, DISP_POSI_CONSS, DISP_STRI_CONSS) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_CUTS);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_CUTS, DISP_DESC_CUTS, DISP_HEAD_CUTS,
            SCIP_DISPSTATUS_AUTO,
            dispCopyDefault,
            NULL, NULL, NULL, NULL, NULL, SCIPdispOutputNAppliedCuts, NULL,
            DISP_WIDT_CUTS, DISP_PRIO_CUTS, DISP_POSI_CUTS, DISP_STRI_CUTS) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_SEPAROUNDS);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_SEPAROUNDS, DISP_DESC_SEPAROUNDS, DISP_HEAD_SEPAROUNDS,
            SCIP_DISPSTATUS_AUTO,
            dispCopyDefault,
            NULL, NULL, NULL, NULL, NULL, SCIPdispOutputNSepaRounds, NULL,
            DISP_WIDT_SEPAROUNDS, DISP_PRIO_SEPAROUNDS, DISP_POSI_SEPAROUNDS, DISP_STRI_SEPAROUNDS) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_POOLSIZE);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_POOLSIZE, DISP_DESC_POOLSIZE, DISP_HEAD_POOLSIZE,
            SCIP_DISPSTATUS_AUTO,
            dispCopyDefault,
            NULL, NULL, NULL, NULL, NULL, SCIPdispOutputCutPoolSize, NULL,
            DISP_WIDT_POOLSIZE, DISP_PRIO_POOLSIZE, DISP_POSI_POOLSIZE, DISP_STRI_POOLSIZE) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_LPOBJ);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_LPOBJ, DISP_DESC_LPOBJ, DISP_HEAD_LPOBJ,
            SCIP_DISPSTATUS_AUTO,
            dispCopyDefault,
            NULL, NULL, NULL, NULL, NULL, SCIPdispOutputLPObjval, NULL,
            DISP_WIDT_LPOBJ, DISP_PRIO_LPOBJ, DISP_POSI_LPOBJ, DISP_STRI_LPOBJ) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_MVARS);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_MVARS, DISP_DESC_MVARS, DISP_HEAD_MVARS,
            SCIP_DISPSTATUS_AUTO, dispCopyDefault, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMvars, NULL,
            DISP_WIDT_MVARS, DISP_PRIO_MVARS, DISP_POSI_MVARS, DISP_STRI_MVARS) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_MCONSS);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_MCONSS, DISP_DESC_MCONSS, DISP_HEAD_MCONSS,
            SCIP_DISPSTATUS_AUTO, dispCopyDefault, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMconss, NULL,
            DISP_WIDT_MCONSS, DISP_PRIO_MCONSS, DISP_POSI_MCONSS, DISP_STRI_MCONSS) );
   }
   tmpdisp = SCIPfindDisp(scip, DISP_NAME_MCUTS);
   if( tmpdisp == NULL )
   {
      SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME_MCUTS, DISP_DESC_MCUTS, DISP_HEAD_MCUTS,
            SCIP_DISPSTATUS_AUTO, dispCopyDefault, NULL, NULL, NULL, NULL, NULL, SCIPdispOutputMcuts, NULL,
            DISP_WIDT_MCUTS, DISP_PRIO_MCUTS, DISP_POSI_MCUTS, DISP_STRI_MCUTS) );
   }
   return SCIP_OKAY;
}
