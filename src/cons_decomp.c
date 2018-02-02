/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

/**@file   cons_decomp.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for structure detection
 * @author Martin Bergner
 * @author Christian Puchert
 *
 * This constraint handler will run all registered structure detectors in
 * increasing priority until the first detector finds a suitable structure.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "cons_decomp.h"
#include "dec_connected.h"
#include "gcg.h"
#include "struct_detector.h"
#include "string.h"
#include "scip_misc.h"
#include "scip/clock.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "decomp"
#define CONSHDLR_DESC          "constraint handler for structure detection"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                          *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

#define DEFAULT_CREATEBASICDECOMP FALSE /**< indicates whether to create a decomposition with all constraints in the master if no other specified */
/*
 * Data structures
 */


/** constraint handler data */
struct SCIP_ConshdlrData
{
   DEC_DECOMP**          decdecomps;         /**< array of decomposition structures */
   DEC_DETECTOR**        detectors;          /**< array of structure detectors */
   int*                  priorities;         /**< priorities of the detectors */
   int                   ndetectors;         /**< number of detectors */
   SCIP_CLOCK*           detectorclock;      /**< clock to measure detection time */
   SCIP_Bool             hasrun;             /**< flag to indicate whether we have already detected */
   int                   ndecomps;           /**< number of decomposition structures  */
   SCIP_Bool             createbasicdecomp;  /**< indicates whether to create a decomposition with all constraints in the master if no other specified */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of constraint handler
 */

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitDecomp)
{ /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->hasrun = FALSE;

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);

      SCIP_CALL( SCIPresetClock(scip, detector->dectime) );
      if( detector->initDetector != NULL )
      {
         SCIPdebugMessage("Calling initDetector of %s\n", detector->name);
         SCIP_CALL( (*detector->initDetector)(scip, detector) );
      }
   }

   return SCIP_OKAY;
}

/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitDecomp)
{ /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(conshdlr != NULL);
   assert(scip != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->ndecomps > 0 )
   {
      for( i = 0; i < conshdlrdata->ndecomps; ++i )
      {
         SCIP_CALL( DECdecompFree(scip, &conshdlrdata->decdecomps[i]) );
      }
      SCIPfreeMemoryArray(scip, &conshdlrdata->decdecomps);
      conshdlrdata->decdecomps = NULL;
      conshdlrdata->ndecomps = 0;
   }
   conshdlrdata->hasrun = FALSE;

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);

      detector->ndecomps = 0;
      SCIPfreeMemoryArrayNull(scip, &detector->decomps);
      if( detector->exitDetector != NULL )
      {
         SCIPdebugMessage("Calling exitDetector of %s\n", detector->name);
         SCIP_CALL( (*detector->exitDetector)(scip, detector) );
      }
   }

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeDecomp)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->detectorclock) );

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);

      if( detector->freeDetector != NULL )
      {
         SCIPdebugMessage("Calling freeDetector of %s\n", detector->name);
         SCIP_CALL( (*detector->freeDetector)(scip, detector) );
      }
      SCIP_CALL( SCIPfreeClock(scip, &detector->dectime) );
      SCIPfreeBlockMemory(scip, &detector);
   }

   /* @todo: This is also done in consExitDecomp() and therefore probably makes no sense here. */
   if( conshdlrdata->ndecomps > 0 )
   {
      for( i = 0; i < conshdlrdata->ndecomps; ++i )
      {
         SCIP_CALL( DECdecompFree(scip, &conshdlrdata->decdecomps[i]) );
      }
      SCIPfreeMemoryArray(scip, &conshdlrdata->decdecomps);
   }

   SCIPfreeMemoryArray(scip, &conshdlrdata->priorities);
   SCIPfreeMemoryArray(scip, &conshdlrdata->detectors);
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->decdecomps);
   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforeDecomp)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpDecomp)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsDecomp)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckDecomp)
{
   /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}



/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockDecomp)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for decomp constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create decomp constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   assert(conshdlrdata != NULL);

   conshdlrdata->decdecomps = NULL;
   conshdlrdata->ndecomps = 0;
   conshdlrdata->ndetectors = 0;
   conshdlrdata->priorities = NULL;
   conshdlrdata->detectors = NULL;
   conshdlrdata->hasrun = FALSE;

   SCIP_CALL( SCIPcreateWallClock(scip, &conshdlrdata->detectorclock) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpDecomp, consEnfopsDecomp, consCheckDecomp, consLockDecomp,
         conshdlrdata) );
   assert(conshdlr != FALSE);

   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforeDecomp) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeDecomp) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitDecomp) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitDecomp) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/decomp/createbasicdecomp", "indicates whether to create a decomposition with all constraints in the master if no other specified", &conshdlrdata->createbasicdecomp, FALSE, DEFAULT_CREATEBASICDECOMP, NULL, NULL) );

   return SCIP_OKAY;
}

/** sets (and adds) the decomposition structure **/
SCIP_RETCODE SCIPconshdlrDecompAddDecdecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->ndecomps == 0 )
   {
      assert(conshdlrdata->decdecomps == NULL);
      SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->decdecomps, 1) ); /*lint !e506*/
      conshdlrdata->decdecomps[0] = decdecomp;
      conshdlrdata->ndecomps = 1;
   }
   else
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->decdecomps, (size_t)conshdlrdata->ndecomps+1) );
      conshdlrdata->decdecomps[conshdlrdata->ndecomps] = conshdlrdata->decdecomps[0];
      conshdlrdata->decdecomps[0] = decdecomp;
      conshdlrdata->ndecomps += 1;
   }
   return SCIP_OKAY;
}


/** returns the decomposition structure **/
DEC_DECOMP** SCIPconshdlrDecompGetDecdecomps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->decdecomps;
}

/** returns the decomposition structure **/
int SCIPconshdlrDecompGetNDecdecomps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->ndecomps;
}

/** returns the data of the provided detector */
DEC_DETECTORDATA* DECdetectorGetData(
   DEC_DETECTOR*         detector            /**< detector data structure */
   )
{
   assert(detector != NULL);
   return detector->decdata;

}

/** returns the name of the provided detector */
const char* DECdetectorGetName(
   DEC_DETECTOR*         detector            /**< detector data structure */
   )
{
   assert(detector != NULL);
   return detector->name;
}

/** searches for the detector and returns it or returns NULL if detector is not found*/
DEC_DETECTOR* DECfindDetector(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of the detector */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
      return NULL;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);
      if( strcmp(detector->name, name) == 0 )
      {
         return detector;
      }
   }

   return NULL;
}

/** includes the detector */
SCIP_RETCODE DECincludeDetector(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the detector */
   const char            decchar,            /**< display character of the detector */
   const char*           description,        /**< description of the detector */
   int                   priority,           /**< priority of the detector */
   SCIP_Bool             enabled,            /**< whether the detector should be enabled by default */
   SCIP_Bool             skip,               /**< whether the detector should be skipped if others found structure */
   DEC_DETECTORDATA*     detectordata,       /**< the associated detector data (or NULL) */
   DEC_DECL_DETECTSTRUCTURE((*detectStructure)), /**< the method that will detect the structure (must not be NULL)*/
   DEC_DECL_FREEDETECTOR((*freeDetector)),   /**< destructor of detector (or NULL) */
   DEC_DECL_INITDETECTOR((*initDetector)),   /**< initialization method of detector (or NULL) */
   DEC_DECL_EXITDETECTOR((*exitDetector))    /**< deinitialization method of detector (or NULL) */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   DEC_DETECTOR *detector;
   char setstr[SCIP_MAXSTRLEN];
   char descstr[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(name != NULL);
   assert(description != NULL);
   assert(detectStructure != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   assert(detectStructure != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &detector) );
   assert(detector != NULL);

   SCIPdebugMessage("Adding detector %i: %s\n", conshdlrdata->ndetectors+1, name);

#ifndef NDEBUG
   assert(DECfindDetector(scip, name) == NULL);
#endif

   detector->decdata = detectordata;
   detector->name = name;
   detector->description = description;
   detector->decchar = decchar;

   detector->freeDetector = freeDetector;
   detector->initDetector = initDetector;
   detector->exitDetector = exitDetector;
   detector->detectStructure = detectStructure;

   detector->priority = priority;
   detector->enabled = enabled;
   detector->skip = skip;
   detector->ndecomps = 0;
   detector->decomps = NULL;
   SCIP_CALL( SCIPcreateWallClock(scip, &(detector->dectime)) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/enabled", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> is enabled", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->enabled), FALSE, enabled, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/skip", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> should be skipped if others found decompositions", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->skip), FALSE, skip, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/priority", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "priority of detector <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->priority), FALSE, priority, INT_MIN, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->detectors, (size_t)conshdlrdata->ndetectors+1) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->priorities,(size_t) conshdlrdata->ndetectors+1) );

   conshdlrdata->detectors[conshdlrdata->ndetectors] = detector;
   conshdlrdata->ndetectors = conshdlrdata->ndetectors+1;

   return SCIP_OKAY;

}

/** returns the remaining time of scip that the decomposition may use */
SCIP_Real DECgetRemainingTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real timelimit;
   assert(scip != NULL);
   SCIP_CALL_ABORT(SCIPgetRealParam(scip, "limits/time", &timelimit));
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   return timelimit;
}

/** interface method to detect the structure */
SCIP_RETCODE DECdetectStructure(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< Result pointer to indicate whether some structure was found */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_Real* scores;
   int i;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT || SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "No problem exists, cannot detect structure!\n");

      if( SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
         conshdlrdata->hasrun = TRUE;

      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
      SCIP_CALL( SCIPtransformProb(scip) );

   SCIP_CALL( SCIPresetClock(scip, conshdlrdata->detectorclock) );
   SCIP_CALL( SCIPstartClock(scip, conshdlrdata->detectorclock) );

   if( conshdlrdata->ndecomps == 0 )
   {
      for( i = 0; i < conshdlrdata->ndetectors; ++i )
      {
         DEC_DETECTOR *detector;
         detector = conshdlrdata->detectors[i];
         assert(detector != NULL);
         conshdlrdata->priorities[i] = detector->priority;
      }

      SCIPdebugMessage("Sorting %i detectors\n", conshdlrdata->ndetectors);
      SCIPsortIntPtr(conshdlrdata->priorities, (void**)conshdlrdata->detectors, conshdlrdata->ndetectors);

      SCIPdebugMessage("Trying %d detectors.\n", conshdlrdata->ndetectors);
      for( i = 0; i < conshdlrdata->ndetectors; ++i )
      {
         DEC_DETECTOR* detector;
         DEC_DECOMP** decdecomps;
         int ndecdecomps;

         ndecdecomps = -1;
         detector = conshdlrdata->detectors[i];
         assert(detector != NULL);
         if( !detector->enabled )
            continue;
         decdecomps = NULL;

         SCIPdebugMessage("Calling detectStructure of %s: ", detector->name);
         SCIP_CALL( SCIPstartClock(scip, detector->dectime) );
         SCIP_CALL( (*detector->detectStructure)(scip, detector->decdata, &decdecomps, &ndecdecomps,  result) );
         SCIP_CALL( SCIPstopClock(scip, detector->dectime) );
         SCIPdebugPrintf("(time %.6f) ", SCIPclockGetTime(detector->dectime));

         if( *result == SCIP_SUCCESS )
         {
            int j;
            assert(ndecdecomps >= 0);
            assert(decdecomps != NULL || ndecdecomps == 0);

            SCIPdebugMessage("We originally have %d decompositions, ", ndecdecomps);
            for( j = 0; j < ndecdecomps; ++j )
            {
               assert(decdecomps != NULL);
               DECdecompSetDetector(decdecomps[j], detector);
            }
            if( ndecdecomps > 2 )
            {
               int nunique = DECfilterSimilarDecompositions(scip, decdecomps, ndecdecomps);

               for( j = nunique; j < ndecdecomps; ++j )
               {
                  SCIP_CALL( DECdecompFree(scip, &(decdecomps[j])) );
                  decdecomps[j] = NULL;
               }

               ndecdecomps = nunique;
            }
            SCIPdebugPrintf("%d after filtering!\n", ndecdecomps);

            SCIP_CALL( SCIPreallocMemoryArray(scip, &(conshdlrdata->decdecomps), ((size_t)conshdlrdata->ndecomps+ndecdecomps)) );
            BMScopyMemoryArray(&(conshdlrdata->decdecomps[conshdlrdata->ndecomps]), decdecomps, ndecdecomps); /*lint !e866*/
            conshdlrdata->ndecomps += ndecdecomps;
            detector->ndecomps = ndecdecomps;
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &detector->decomps, decdecomps, detector->ndecomps) );
         }
         else
         {
            SCIPdebugPrintf("Failure!\n");
         }
         SCIPfreeMemoryArrayNull(scip, &decdecomps);
      }
   }
   else
   {
      SCIP_CALL( DECdecompTransform(scip, conshdlrdata->decdecomps[0]) );
   }
   /* evaluate all decompositions and sort them by score */
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, conshdlrdata->ndecomps) );
   for( i = 0; i < conshdlrdata->ndecomps; ++i )
   {
      DEC_SCORES score;
      score.totalscore = 0.0;

      SCIP_CALL( DECevaluateDecomposition(scip, conshdlrdata->decdecomps[i], &score) );
      scores[i] = score.totalscore;
   }

   SCIPsortRealPtr(scores, (void**)conshdlrdata->decdecomps, conshdlrdata->ndecomps);
   SCIPfreeBufferArray(scip, &scores);

   SCIP_CALL( SCIPstopClock(scip, conshdlrdata->detectorclock) );

   if( conshdlrdata->ndecomps > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Chosen decomposition with %d blocks of type %s.\n",
         DECdecompGetNBlocks(conshdlrdata->decdecomps[0]), DECgetStrType(DECdecompGetType(conshdlrdata->decdecomps[0])));
      GCGsetStructDecdecomp(scip, conshdlrdata->decdecomps[0]);
      *result = SCIP_SUCCESS;
   }
   else
   {
      *result = SCIP_DIDNOTFIND;
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "No decomposition found!\n");
   }
   SCIPdebugMessage("Detection took %fs\n", SCIPclockGetTime(conshdlrdata->detectorclock));

   /* show that we done our duty */
   conshdlrdata->hasrun = TRUE;

   return SCIP_OKAY;
}

/** write out all detected or provided decompositions */
SCIP_RETCODE DECwriteAllDecomps(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 directory,          /**< directory for decompositions */
   char*                 extension           /**< extension for decompositions */
   )
{
   char name[SCIP_MAXSTRLEN];
   char outname[SCIP_MAXSTRLEN];
   char *pname;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   DEC_DETECTOR *detector;
   DEC_DECOMP *decomp;
   DEC_DECOMP *tmp;
   int i;
   int j;

   assert(scip != NULL);
   assert(extension != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->ndecomps == 0 )
   {
      SCIPwarningMessage(scip, "No decomposition available.\n");
      return SCIP_OKAY;
   }

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   SCIPsplitFilename(name, NULL, &pname, NULL, NULL);

   tmp = conshdlrdata->decdecomps[0];

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      detector =  conshdlrdata->detectors[i];
      assert(detector != NULL);

      for( j = 0; j < detector->ndecomps; ++j )
      {
         decomp = detector->decomps[j];
         assert(decomp != NULL);
         if( directory != NULL )
         {
            (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s/%s_%c_%d_%d.%s", directory, pname, detector->decchar, DECdecompGetNBlocks(decomp), j, extension);
         }
         else
         {
            (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%c_%d_%d.%s", pname, detector->decchar, DECdecompGetNBlocks(decomp), j, extension);

         }
         conshdlrdata->decdecomps[0] = decomp;
         SCIP_CALL( SCIPwriteTransProblem(scip, outname, extension, FALSE) );
      }
   }

   /** further, get all read in decompositions */
   for( i = 0; i < conshdlrdata->ndecomps; ++i )
   {
      decomp = conshdlrdata->decdecomps[i];
      detector =  DECdecompGetDetector(decomp);

      if( detector != NULL )
         continue;

      (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%d.%s", pname, DECdecompGetNBlocks(decomp), extension);

      conshdlrdata->decdecomps[0] = decomp;
      SCIP_CALL( SCIPwriteTransProblem(scip, outname, extension, FALSE) );
   }

   conshdlrdata->decdecomps[0] = tmp;

   return SCIP_OKAY;
}

/** returns the best known decomposition, if available and NULL otherwise */
DEC_DECOMP* DECgetBestDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->ndecomps > 0 )
      return conshdlrdata->decdecomps[0];

   else if ( conshdlrdata->createbasicdecomp)
   {
      SCIP_RETCODE retcode;
      DEC_DECOMP* decomp = NULL;
      retcode = DECcreateBasicDecomp(scip, &decomp);
      assert(retcode == SCIP_OKAY);
      assert(decomp != NULL );

      retcode = SCIPconshdlrDecompAddDecdecomp(scip, decomp);
      if( retcode != SCIP_OKAY )
      {
         SCIPerrorMessage("Could not add decomp to cons_decomp!\n");
         return NULL;
      }

      assert(conshdlrdata->ndecomps > 0);
      assert(conshdlrdata->decdecomps[0] != NULL);
      return conshdlrdata->decdecomps[0];
   }

   return NULL;
}

/** writes out a list of all detectors */
void DECprintListOfDetectors(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int ndetectors;
   int i;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   ndetectors = conshdlrdata->ndetectors;

   SCIPdialogMessage(scip, NULL, " detector             priority char  description\n --------------       -------- ----  -----------\n");

   for( i = 0; i < ndetectors; ++i )
   {
      SCIPdialogMessage(scip, NULL,  " %-20s", conshdlrdata->detectors[i]->name);
      SCIPdialogMessage(scip, NULL,  " %8d    %c ", conshdlrdata->detectors[i]->priority, conshdlrdata->detectors[i]->decchar);
      SCIPdialogMessage(scip, NULL,  " %s\n", conshdlrdata->detectors[i]->description);
   }
}

/** returns whether the detection has been performed */
SCIP_Bool DEChasDetectionRun(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->hasrun;
}

/** returns the character of the detector */
char DECdetectorGetChar(
   DEC_DETECTOR*         detector            /**< pointer to detector */
)
{
   if( detector == NULL )
     return '0';
   else
      return detector->decchar;
}

/** display statistics about detectors */
SCIP_RETCODE GCGprintDetectorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   int j;
   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(scip != NULL);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Detector statistics:       time     number     blocks\n");
   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  %-10.10s       :   %8.2f %10d    ", conshdlrdata->detectors[i]->name, SCIPclockGetTime(conshdlrdata->detectors[i]->dectime), conshdlrdata->detectors[i]->ndecomps );
      for( j = 0; j < conshdlrdata->detectors[i]->ndecomps; ++j )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " %d", DECdecompGetNBlocks(conshdlrdata->detectors[i]->decomps[j]));
      }
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "\n");
   }
   return SCIP_OKAY;
}

/** resets the parameters to their default value */
static
SCIP_RETCODE setDetectionDefault(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{ /*lint --e{715}*/
   int i;
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      char paramname[SCIP_MAXSTRLEN];
      SCIP_Bool paramval;
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/enabled", conshdlrdata->detectors[i]->name);

      SCIP_CALL( SCIPresetParam(scip, paramname) );
      if( !quiet )
      {
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
      }
   }

   return SCIP_OKAY;
}

/** sets the parameters to aggressive values */
static
SCIP_RETCODE setDetectionAggressive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{ /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   return SCIP_OKAY;
}

/** disables detectors */
static
SCIP_RETCODE setDetectionOff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{ /*lint --e{715}*/
   int i;
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      char paramname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/enabled", conshdlrdata->detectors[i]->name);

      SCIP_CALL( SCIPsetBoolParam(scip, paramname, FALSE) );
      if( !quiet )
      {
         SCIPinfoMessage(scip, NULL, "%s = FALSE\n", paramname);
      }
   }

   return SCIP_OKAY;
}

/** sets the parameters to fast values */
static
SCIP_RETCODE setDetectionFast(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{ /*lint --e{715} */
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   return SCIP_OKAY;
}

/** sets detector parameters values to
 *
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all detector parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for detection is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the detectors produce more decompositions
 *  - SCIP_PARAMSETTING_OFF which turns off all detection
 */
SCIP_RETCODE GCGsetDetection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(scip != NULL);

   switch( paramsetting )
   {
   case SCIP_PARAMSETTING_AGGRESSIVE:
      SCIP_CALL( setDetectionDefault(scip, conshdlrdata, quiet) );
      SCIP_CALL( setDetectionAggressive(scip, conshdlrdata, quiet) );
      break;
   case SCIP_PARAMSETTING_OFF:
      SCIP_CALL( setDetectionOff(scip, conshdlrdata, quiet) );
      break;
   case SCIP_PARAMSETTING_FAST:
      SCIP_CALL( setDetectionDefault(scip, conshdlrdata, quiet) );
      SCIP_CALL( setDetectionFast(scip, conshdlrdata, quiet) );
      break;
   case SCIP_PARAMSETTING_DEFAULT:
      SCIP_CALL( setDetectionDefault(scip, conshdlrdata, quiet) );
      break;
   default:
      SCIPerrorMessage("The given paramsetting is invalid!\n");
      break;
   }

   return SCIP_OKAY;
}
