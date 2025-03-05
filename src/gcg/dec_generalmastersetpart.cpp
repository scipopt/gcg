/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file   dec_generalmastersetpart.cpp
 * 
 * @brief  detector for set partitioning constraints
 * @author Martin Bergner
 *
 * This detector sets the following constraints to master:
 * - set partitioning constraints
 * - constraints with the same nonnegative rhs and lhs
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/dec_generalmastersetpart.h"
#include "gcg/cons_decomp.h"
#include "gcg/class_partialdecomp.h"
#include "gcg/class_detprobdata.h"
#include "gcg/gcg.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"
#include "gcg/scip_misc.h"
#include "scip/clock.h"

#include <iostream>

/* constraint handler properties */
#define DEC_NAME                  "generalmastersetpart"       /**< name of detector */
#define DEC_DESC                  "detector generalmastersetpart" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          0     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  0     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0           /**< priority of the constraint handler for separation */
#define DEC_DECCHAR               '?'         /**< display character of detector */
#define DEC_ENABLED               TRUE        /**< should the detection be enabled */
#define DEC_ENABLEDFINISHING      FALSE       /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the postprocessing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if other detectors found decompositions */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated partialdec */


/*
 * Data structures
 */

/** @todo fill in the necessary detector data */

/** detector handler data */
struct GCG_DetectorData
{
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
#define freeGeneralmastersetpart NULL

/** destructor of detector to free detector data (called before the solving process begins) */
#define exitGeneralmastersetpart NULL

/** detection initialization function of detector (called before solving is about to begin) */
#define initGeneralmastersetpart NULL


static GCG_DECL_PROPAGATEPARTIALDEC(propagatePartialdecGeneralmastersetpart)
{
   *result = SCIP_DIDNOTFIND;
   char decinfo[SCIP_MAXSTRLEN];
   SCIP* origprob = GCGgetOrigprob(gcg);
   SCIP_CLOCK* temporaryClock;
   SCIP_CALL_ABORT(SCIPcreateClock(origprob, &temporaryClock) );
   SCIP_CALL_ABORT( SCIPstartClock(origprob, temporaryClock) );

   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int nvars;
   bool relevant = true;
   SCIP_Real value;

   gcg::PARTIALDECOMP* partialdec = partialdecdetectiondata->workonpartialdec;
   auto& openconss = partialdec->getOpenconss();
   for( auto itr = openconss.cbegin(); itr != openconss.cend(); )
   {
      bool found = false;
      cons = partialdecdetectiondata->detprobdata->getCons(*itr);
      /* set open setpartitioning constraints to master */
      if( GCGconsGetType(origprob, cons) == setpartitioning )
      {
         itr = partialdec->fixConsToMaster(itr);
         found = true;
      }
      /* set constraints with the same nonnegative lhs and rhs to master */
      else if( GCGconsGetType(origprob, cons) != logicor && GCGconsGetType(origprob, cons) != setcovering && GCGconsGetType(origprob, cons) != setpacking )
      {
         nvars = GCGconsGetNVars(origprob, cons);
         vars = NULL;
         vals = NULL;
         value = GCGconsGetLhs(origprob, cons);

         if( SCIPisNegative(origprob, value) || !SCIPisFeasEQ(origprob, GCGconsGetRhs(origprob,cons), value) )
            relevant = false;
         if( nvars > 0 )
         {
            SCIP_CALL( SCIPallocMemoryArray(origprob, &vars, nvars) );
            SCIP_CALL( SCIPallocMemoryArray(origprob, &vals, nvars) );
            SCIP_CALL( GCGconsGetVars(origprob, cons, vars, nvars) );
            SCIP_CALL( GCGconsGetVals(origprob, cons, vals, nvars) );
         }
         for( int j = 0; j < nvars && relevant; ++j )
         {
            assert(vars != NULL);
            assert(vals != NULL);

            if( !SCIPvarIsIntegral(vars[j]) && !SCIPvarIsBinary(vars[j]) )
            {
               SCIPdebugPrintf("(%s is not integral) ", SCIPvarGetName(vars[j]) );
               relevant = false;
            }
            if( !SCIPisEQ(origprob, vals[j], 1.0) )
            {
               SCIPdebugPrintf("(coeff for var %s is %.2f != 1.0) ", SCIPvarGetName(vars[j]), vals[j] );
               relevant = false;
            }
         }
         SCIPfreeMemoryArrayNull(origprob, &vals);
         SCIPfreeMemoryArrayNull(origprob, &vars);

         if( relevant )
         {
            itr = partialdec->fixConsToMaster(itr);
            found = true;
         }
      }
      if( !found )
      {
         ++itr;
      }
   }

   partialdec->sort();
   SCIP_CALL_ABORT( SCIPstopClock(origprob, temporaryClock) );

   partialdecdetectiondata->detectiontime = SCIPgetClockTime(origprob, temporaryClock);
   SCIP_CALL( SCIPallocMemoryArray(origprob, &(partialdecdetectiondata->newpartialdecs), 1) );
   partialdecdetectiondata->newpartialdecs[0] = partialdec;
   partialdecdetectiondata->nnewpartialdecs = 1;
   (void) SCIPsnprintf(decinfo, SCIP_MAXSTRLEN, "genmastersetpart");
   partialdecdetectiondata->newpartialdecs[0]->addDetectorChainInfo(decinfo);
   partialdecdetectiondata->newpartialdecs[0]->addClockTime(SCIPgetClockTime(origprob, temporaryClock));
   // we used the provided partialdec -> prevent deletion
   partialdecdetectiondata->workonpartialdec = NULL;
   SCIP_CALL_ABORT(SCIPfreeClock(origprob, &temporaryClock) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

#define finishPartialdecGeneralmastersetpart NULL
#define detectorPostprocessPartialdecGeneralmastersetpart NULL

static
GCG_DECL_SETPARAMAGGRESSIVE(setParamAggressiveGeneralmastersetpart)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   int newval;
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, TRUE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE ) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxcallround", name);
   SCIP_CALL( SCIPgetIntParam(origprob, setstr, &newval) );
   ++newval;
   SCIP_CALL( SCIPsetIntParam(origprob, setstr, newval ) );
   SCIPinfoMessage(origprob, NULL, "%s = %d\n", setstr, newval);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origmaxcallround", name);
   SCIP_CALL( SCIPgetIntParam(origprob, setstr, &newval) );
   ++newval;
   SCIP_CALL( SCIPsetIntParam(origprob, setstr, newval ) );
   SCIPinfoMessage(origprob, NULL, "%s = %d\n", setstr, newval);

   return SCIP_OKAY;
}


static
GCG_DECL_SETPARAMDEFAULT(setParamDefaultGeneralmastersetpart)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, DEC_ENABLED) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, DEC_ENABLEDFINISHING ) );

   return SCIP_OKAY;
}

static
GCG_DECL_SETPARAMFAST(setParamFastGeneralmastersetpart)
{
   char setstr[SCIP_MAXSTRLEN];
   const char* name = GCGdetectorGetName(detector);
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   SCIP_CALL( SCIPsetBoolParam(origprob, setstr, FALSE ) );

   return SCIP_OKAY;
}



/*
 * detector specific interface methods
 */

/** creates the handler for generalmastersetpart detector and includes it in SCIP */
SCIP_RETCODE GCGincludeDetectorGeneralmastersetpart(
   GCG*                 gcg                  /**< GCG data structure */
   )
{
   GCG_DETECTORDATA* detectordata;

   /**@todo create generalmastersetpart detector data here*/
   detectordata = NULL;

   SCIP_CALL(
      GCGincludeDetector(gcg, DEC_NAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, detectordata, freeGeneralmastersetpart, initGeneralmastersetpart, exitGeneralmastersetpart, propagatePartialdecGeneralmastersetpart, finishPartialdecGeneralmastersetpart, detectorPostprocessPartialdecGeneralmastersetpart, setParamAggressiveGeneralmastersetpart, setParamDefaultGeneralmastersetpart, setParamFastGeneralmastersetpart));

   /**@todo add generalmastersetpart detector parameters */

   return SCIP_OKAY;
}





