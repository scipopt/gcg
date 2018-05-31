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

/**@file   dec_random.c
 * @ingroup DETECTORS
 * @brief  Random structure detection for testing purposes
 * @author Martin Bergner
 *
 * This detector will partition the constraints of the problem randomly.
 * For each constraint, it will randomly pick a number between 0 and the
 * maxblocks parameter. Constraints assigned to maxblocks will be put in the
 * master problem. You can set the random seed via the seed parameter.
 *
 * If that maxblocks parameter is set to -1, it will default to number of
 * constraints.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dec_random.h"
#include "cons_decomp.h"
#include "pub_decomp.h"
#include "scip/clock.h"
#include "scip/misc.h"
#include <string.h>

/* constraint handler properties */
#define DEC_DETECTORNAME          "random"    /**< name of detector */
#define DEC_DESC                  "Random structure detection" /**< description of detector */
#define DEC_PRIORITY              -10         /**< priority of the constraint handler for separation */
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /** last round the detector gets called                              */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                              */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_DECCHAR               'R'         /**< display character of detector */
#define DEC_ENABLED               FALSE       /**< should the detection be enabled */
#define DEC_ENABLEDORIGINAL       FALSE       /**< should the detection of the original problem be enabled */
#define DEC_ENABLEDFINISHING      FALSE       /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the postprocessing be enabled */
#define DEC_SKIP                  FALSE       /**< should detector be skipped if others found detections */
#define DEC_USEFULRECALL          FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */
#define DEC_LEGACYMODE            FALSE       /**< should (old) DETECTSTRUCTURE method also be used for detection */

#define DEFAULT_MAXBLOCKS        -1          /**< the maximal number of blocks, -1 defaults to average number of constraints */
#define DEFAULT_AVGCONSPERBLOCK  100         /**< average constraints per block to limit the maximal block number */
#define DEFAULT_RANDSEED         23          /**< initial random seed */

/*
 * Data structures
 */

/** @todo fill in the necessary detector data */

/** detector handler data */
struct DEC_DetectorData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   int                   maxblocks;          /**< the maximal number of blocks, -1 defaults to nconss/maxconsperblock */
   int                   avgconsperblock;    /**< the average number of constraints per block */
   SCIP_HASHMAP*         constoblock;        /**< hashmap to store partition */
   int                   nblocks;            /**< number of actual blocks found */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** randomly assings constraints to blocks or the master */
static
SCIP_RETCODE findRandomPartition(
   SCIP*              scip,                  /**< SCIP data structure */
   DEC_DETECTORDATA*  detectordata           /**< detector data structure */
   )
{
   int nconss;
   SCIP_CONS** conss;
   int i;
   int maxblocks;
   int* consblocks;
   int oldblock;
   int nblocks;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(detectordata->constoblock != NULL);

   nconss = SCIPgetNConss(scip);
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &conss, SCIPgetConss(scip), nconss) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &consblocks, nconss) );
   BMSclearMemoryArray(consblocks, nconss);
   
   if( detectordata->maxblocks == -1 )
      maxblocks = nconss/detectordata->avgconsperblock;
   else
      maxblocks = detectordata->maxblocks;

   nblocks = 0;

   for( i = 0; i < nconss; ++i )
      consblocks[i] = SCIPrandomGetInt(detectordata->randnumgen, 0, maxblocks);

   SCIPsortIntPtr(consblocks, (void**)conss, nconss);
   oldblock = -1;
   for( i = 0; i < nconss; ++i )
   {
      assert(consblocks[i] >= oldblock);
      if( consblocks[i] != oldblock )
      {
         oldblock = consblocks[i];
         ++nblocks;
      }
      SCIPdebugMessage("Assigning cons <%s> to block %d.\n", SCIPconsGetName(conss[i]), nblocks);
      assert(nblocks > 0);
      SCIP_CALL( SCIPhashmapInsert(detectordata->constoblock, conss[i], (void*)(size_t) (nblocks)) ); /*lint !e866*/
   }

   detectordata->nblocks = nblocks;

   SCIPfreeMemoryArray(scip, &consblocks);
   SCIPfreeMemoryArray(scip, &conss);
   return SCIP_OKAY;
}


/*
 * detector callback methods
 */

/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(detectorFreeRandom)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** detector initialization method (called after problem was transformed) */
static
DEC_DECL_INITDETECTOR(detectorInitRandom)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   detectordata->maxblocks = 0;
   detectordata->constoblock = NULL;
   detectordata->nblocks = 0;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &detectordata->randnumgen,
         SCIPinitializeRandomSeed(scip, DEFAULT_RANDSEED), TRUE) );

   return SCIP_OKAY;
}

/** detector deinitialization method (called before the transformed problem is freed) */
static
DEC_DECL_EXITDETECTOR(detectorExitRandom)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   /* free random number generator */

   SCIPfreeRandom(scip, &detectordata->randnumgen);

   return SCIP_OKAY;
}

/** detection function of detector */
static
DEC_DECL_DETECTSTRUCTURE(detectorDetectRandom)
{ /*lint --e{715}*/
   *result = SCIP_DIDNOTFIND;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting random structure:");

   detectordata->nblocks = 0;

   SCIP_CALL( SCIPhashmapCreate(&detectordata->constoblock, SCIPblkmem(scip), SCIPgetNConss(scip)) );

   SCIP_CALL( findRandomPartition(scip, detectordata) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " found %d blocks.\n", detectordata->nblocks);

   if( detectordata->nblocks > 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, 1) ); /*lint !e506*/
      SCIP_CALL( DECdecompCreate(scip, &((*decdecomps)[0])) );

      SCIP_CALL( DECfilloutDecompFromConstoblock(scip, (*decdecomps)[0], detectordata->constoblock, detectordata->nblocks, FALSE) );

      /* delete, debugging */
      SCIP_CALL( DECdecompCheckConsistency( scip, (*decdecomps)[0] ) );

      *ndecdecomps = 1;

      *result = SCIP_SUCCESS;
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " not found.\n");
   }

   /* do not free hashmap since this would also delete assignments in decdecomps */
   // SCIPhashmapFree(&detectordata->constoblock);
   detectordata->constoblock = NULL;

   return SCIP_OKAY;
}

#define detectorPropagateSeeedRandom NULL
#define detectorFinishSeeedRandom NULL
#define detectorPostprocessSeeedRandom NULL


#define setParamAggressiveRandom NULL
#define setParamDefaultRandom NULL
#define setParamFastRandom NULL

/*
 * detector specific interface methods
 */

/** creates the handler for random detector and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectorRandom(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   detectordata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLEDORIGINAL, DEC_ENABLEDFINISHING,DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, DEC_LEGACYMODE,
      detectordata, detectorDetectRandom, detectorFreeRandom, detectorInitRandom, detectorExitRandom, detectorPropagateSeeedRandom, NULL, NULL, detectorFinishSeeedRandom, detectorPostprocessSeeedRandom, setParamAggressiveRandom, setParamDefaultRandom, setParamFastRandom) );


   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/random/maxblocks", "the maximal number of blocks, -1 defaults to avgconsperblock",
      &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, -1, INT_MAX, NULL, NULL ) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/detectors/random/avgconsperblock", "average constraints per block",
      &detectordata->avgconsperblock, FALSE, DEFAULT_AVGCONSPERBLOCK, 1, 10000, NULL, NULL ) );

   return SCIP_OKAY;
}
